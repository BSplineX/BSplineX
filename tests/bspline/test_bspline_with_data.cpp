#include <fstream>

// Third-party includes
#include "../../cmake/test_cmake_inclusion/third_party/bsplinex-src/include/BSplineX/bspline/bspline_types.hpp"
#include "../../third_party/catch2-src/extras/catch_amalgamated.hpp"
#include "BSplineX/bspline/bspline.hpp"
#include "matchers.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>
#include <nlohmann/json.hpp>

using json   = nlohmann::json;
using real_t = double;

// BSplineX includes

using namespace Catch::Matchers;
using namespace bsplinex;
using namespace bsplinex::bspline;

TEST_CASE("test", "[test]")
{
  std::ifstream f("/home/gabriele/Documents/BSplineX/BSplineX/tests/bspline/data/open.json");
  json data = json::parse(f);

  size_t index  = 4;
  size_t degree = data[index]["bspline"]["degree"].get<size_t>();
  types::OpenNonUniform<real_t> bspline{
      data[index]["bspline"]["knots"].get<std::vector<real_t>>(),
      data[index]["bspline"]["ctrl"].get<std::vector<real_t>>(),
      degree,
  };

  auto x_eval = data[index]["x_eval"].get<std::vector<real_t>>();

  SECTION("evaluate(...)")
  {
    auto y_eval = data[index]["bspline"]["y_eval"].get<std::vector<real_t>>();
    for (size_t i{index}; i < x_eval.size(); i++)
    {
      REQUIRE_THAT(bspline.evaluate(x_eval.at(i)), WithinRel(y_eval.at(i)));
    }
    // test also vectorized evaluation
    REQUIRE_THAT(bspline.evaluate(x_eval), VectorsWithinRel(y_eval));
    REQUIRE_THAT(bspline.evaluate(x_eval), VectorsWithinAbsRel(y_eval));
    REQUIRE_THROWS_WITH(
        bspline.evaluate(*std::min_element(x_eval.begin(), x_eval.end()) - 1),
        "Extrapolation explicitly set to NONE"
    );
  }

  SECTION("fit(...)")
  {
    auto [domain_left, domain_right] =
        data[index]["bspline"]["domain"].get<std::pair<real_t, real_t>>();
    std::vector<real_t> x = data[index]["x"].get<std::vector<real_t>>();
    std::vector<real_t> y = data[index]["y"].get<std::vector<real_t>>();
    std::vector<real_t> x_fit, y_fit;
    x_fit.reserve(x.size());
    y_fit.reserve(y.size());
    for (size_t i{0}; i < x.size(); i++)
    {
      if (domain_left <= x.at(i) && x.at(i) <= domain_right)
      {
        x_fit.push_back(x.at(i));
        y_fit.push_back(y.at(i));
      }
    }
    auto y_eval_fit   = data[index]["bspline_fit"]["y_eval"].get<std::vector<real_t>>();
    auto knots_fit    = data[index]["bspline_fit"]["knots"].get<std::vector<real_t>>();
    auto ctrl_pts_fit = data[index]["bspline_fit"]["ctrl"].get<std::vector<real_t>>();

    bspline.fit(x_fit, y_fit);

    REQUIRE_THAT(bspline.get_knots(), VectorsWithinAbsRel(knots_fit));
    REQUIRE_THAT(bspline.get_control_points(), VectorsWithinAbsRel(ctrl_pts_fit));
    REQUIRE_THAT(bspline.evaluate(x_eval), VectorsWithinAbsRel(y_eval_fit));
  }
}
