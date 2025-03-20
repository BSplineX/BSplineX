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

  types::OpenNonUniform<real_t> bspline{
      data[0]["bspline"]["knots"].get<std::vector<real_t>>(),
      data[0]["bspline"]["ctrl"].get<std::vector<real_t>>(),
      data[0]["bspline"]["degree"].get<size_t>(),
  };

  auto x_eval = data[0]["x_eval"].get<std::vector<real_t>>();

  SECTION("evaluate(...)")
  {
    auto y_eval = data[0]["bspline"]["y_eval"].get<std::vector<real_t>>();
    for (size_t i{0}; i < x_eval.size(); i++)
    {
      REQUIRE_THAT(bspline.evaluate(x_eval.at(i)), WithinRel(y_eval.at(i)));
    }
    // test also vectorized evaluation
    REQUIRE_THAT(bspline.evaluate(x_eval), VectorsWithinRel(y_eval));
  }
}
