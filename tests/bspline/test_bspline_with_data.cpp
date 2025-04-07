// Standard includes
#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

// Third-party includes
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <nlohmann/json.hpp>
#include <nlohmann/json_fwd.hpp>

// BSplineX includes
#include "BSplineX/bspline/bspline_lsq.hpp"
#include "BSplineX/bspline/bspline_types.hpp"
#include "BSplineX/types.hpp"
#include "matchers.hpp"

using real_t = double;

using namespace Catch::Matchers;
using namespace bsplinex;
using namespace bsplinex::bspline;

// clang-format off
#define BSPLINE_TEST_TYPES \
types::OpenUniform<real_t>, \
types::OpenNonUniform<real_t>, \
types::ClampedUniform<real_t>, \
types::ClampedNonUniform<real_t>, \
types::PeriodicUniform<real_t>, \
types::PeriodicNonUniform<real_t>
// clang-format on

namespace
{
nlohmann::json load_test_data(std::filesystem::path const &filepath)
{
  std::ifstream file(filepath);
  REQUIRE(file.is_open());
  return nlohmann::json::parse(file);
}

std::string get_test_name(
    std::string const &boundary_condition,
    std::string const &curve_type,
    size_t degree,
    size_t num_knots
)
{
  return "type=" + boundary_condition + "_" + curve_type + "_degree=" + std::to_string(degree) +
         "_knots=" + std::to_string(num_knots);
}

template <typename BSplineType>
std::string get_data_path()
{
  std::filesystem::path path(TEST_RESOURCES_DIR);
  if constexpr (BoundaryCondition::OPEN == BSplineType::boundary_condition_type)
  {
    path /= "open";
  }
  if constexpr (BoundaryCondition::CLAMPED == BSplineType::boundary_condition_type)
  {
    path /= "clamped";
  }
  if constexpr (BoundaryCondition::PERIODIC == BSplineType::boundary_condition_type)
  {
    path /= "periodic";
  }
  // curve
  if constexpr (Curve::UNIFORM == BSplineType::curve_type)
  {
    path /= "uniform";
  }
  if constexpr (Curve::NON_UNIFORM == BSplineType::curve_type)
  {
    path /= "non-uniform";
  }
  path += ".json";
  return path;
}

template <typename BSplineType>
BSplineType build_bspline(std::vector<real_t> knots, std::vector<real_t> ctrl_pts, size_t degree)
{
  if (BoundaryCondition::PERIODIC == BSplineType::boundary_condition_type)
  {
    ctrl_pts.resize(ctrl_pts.size() - degree);
  }
  if constexpr (Curve::UNIFORM == BSplineType::curve_type)
  {
    real_t knots_begin = knots.front();
    real_t knots_end   = knots.back();
    size_t num_knots   = knots.size();

    if (BoundaryCondition::CLAMPED == BSplineType::boundary_condition_type ||
        BoundaryCondition::PERIODIC == BSplineType::boundary_condition_type)
    {
      knots_begin = knots[degree];
      knots_end   = knots[knots.size() - degree - 1];
      num_knots   = knots.size() - 2 * degree;
    }

    return BSplineType({knots_begin, knots_end, num_knots}, {ctrl_pts}, degree);
  }
  else if constexpr (Curve::NON_UNIFORM == BSplineType::curve_type)
  {
    if (BoundaryCondition::CLAMPED == BSplineType::boundary_condition_type ||
        BoundaryCondition::PERIODIC == BSplineType::boundary_condition_type)
    {
      knots = std::vector(std::next(knots.begin(), degree), std::prev(knots.end(), degree));
    }
    return BSplineType({knots}, {ctrl_pts}, degree);
  }
  else
  {
    static_assert(false);
  }
}
} // namespace

TEMPLATE_TEST_CASE("BSpline", "[bspline][template][product]", BSPLINE_TEST_TYPES)
{
  using BSplineType                     = TestType;
  std::filesystem::path const data_path = get_data_path<BSplineType>();
  nlohmann::json const data             = load_test_data(data_path);
  for (auto const &test_data : data)
  {
    auto const boundary_condition = test_data["bspline"]["boundary_condition"].get<std::string>();
    auto const curve_type         = test_data["bspline"]["curve"].get<std::string>();
    auto const degree             = test_data["bspline"]["degree"].get<size_t>();
    auto const num_knots          = test_data["bspline"]["knots"].size();
    auto const test_name = get_test_name(boundary_condition, curve_type, degree, num_knots);

    SECTION(test_name)
    {
      auto bspline = build_bspline<BSplineType>(
          {test_data["bspline"]["knots"].get<std::vector<real_t>>()},
          {test_data["bspline"]["ctrl"].get<std::vector<real_t>>()},
          degree
      );

      auto x_eval = test_data["x_eval"].get<std::vector<real_t>>();
      SECTION("evaluate(...)")
      {
        SECTION("evaluate(..., derivative_order=0)")
        {
          auto y_eval = test_data["bspline"]["y_eval"].get<std::vector<real_t>>();
          for (size_t i{0}; i < x_eval.size(); i++)
          {
            auto matcher = WithinAbsRel(y_eval.at(i));
            REQUIRE_THAT(bspline.evaluate(x_eval.at(i)), matcher);
          }
          // test also vectorized evaluation
          REQUIRE_THAT(bspline.evaluate(x_eval), VectorsWithinAbsRel(y_eval));
          if constexpr (BoundaryCondition::PERIODIC != BSplineType::boundary_condition_type)
          {
            REQUIRE_THROWS_WITH(
                bspline.evaluate(*std::min_element(x_eval.begin(), x_eval.end()) - 1),
                "Extrapolation explicitly set to NONE"
            );
          }
        }
        for (size_t derivative_order{1}; derivative_order <= test_data["derivatives"].size();
             derivative_order++)
        {
          SECTION("evaluate(..., derivative_order=" + std::to_string(derivative_order) + ")")
          {
            auto y_eval =
                test_data["derivatives"][derivative_order - 1]["y_eval"].get<std::vector<real_t>>();
            REQUIRE_THAT(bspline.evaluate(x_eval, derivative_order), VectorsWithinAbsRel(y_eval));
            for (size_t i{0}; i < x_eval.size(); i++)
            {
              auto matcher = WithinAbsRel(y_eval.at(i));
              REQUIRE_THAT(bspline.evaluate(x_eval.at(i), derivative_order), matcher);
            }
            if constexpr (BoundaryCondition::PERIODIC != BSplineType::boundary_condition_type)
            {
              REQUIRE_THROWS_WITH(
                  bspline.evaluate(
                      *std::min_element(x_eval.begin(), x_eval.end()) - 1, derivative_order
                  ),
                  "Extrapolation explicitly set to NONE"
              );
            }
          }
        }
        SECTION("evaluate(..., invalid derivative_order)")
        {
          REQUIRE_THROWS_WITH(
              bspline.evaluate(x_eval[0], degree + 1), "derivative_order must be in [1, degree]"
          );
          REQUIRE_THROWS_WITH(
              bspline.evaluate(x_eval, degree + 1), "derivative_order must be in [1, degree]"
          );
        }
      }

      SECTION("derivative(...)")
      {
        for (size_t derivative_order{1}; derivative_order <= test_data["derivatives"].size();
             derivative_order++)
        {
          SECTION("derivative(derivative_order=" + std::to_string(derivative_order) + ")")
          {
            auto derivative = bspline.derivative(derivative_order);
            REQUIRE(
                derivative.get_degree() ==
                test_data["derivatives"][derivative_order - 1]["degree"].get<size_t>()
            );
            auto y_eval =
                test_data["derivatives"][derivative_order - 1]["y_eval"].get<std::vector<real_t>>();
            REQUIRE_THAT(derivative.evaluate(x_eval), VectorsWithinAbsRel(y_eval));
          }
        }
        SECTION("derivative(invalid derivative_order)")
        {
          REQUIRE_THROWS_WITH(
              bspline.derivative(degree + 1), "derivative_order must be in [1, degree]"
          );
        }
      }

      SECTION("nnz_basis(...)")
      {
        std::vector<real_t> nnz_basis(degree + 1);
        for (size_t derivative_order = 0;
             derivative_order < test_data["bspline"]["nnz_basis"].size();
             ++derivative_order)
        {
          SECTION("nnz_basis(..., derivative_order=" + std::to_string(derivative_order) + ")")
          {
            for (size_t i = 0; i < x_eval.size(); ++i)
            {
              auto [ref_index, ref_nnz_basis] =
                  test_data["bspline"]["nnz_basis"][derivative_order][i]
                      .get<std::pair<size_t, std::vector<real_t>>>();
              std::fill(nnz_basis.begin(), nnz_basis.end(), 0.0);
              auto [index, nnz_basis] = bspline.nnz_basis(x_eval.at(i), derivative_order);
              REQUIRE(index == ref_index);
              REQUIRE_THAT(nnz_basis, VectorsWithinAbsRel(ref_nnz_basis));
            }
          }
        }
      }

      SECTION("fit(...)")
      {
        auto [domain_left, domain_right] =
            test_data["bspline"]["domain"].get<std::pair<real_t, real_t>>();
        std::vector<real_t> x = test_data["x"].get<std::vector<real_t>>();
        std::vector<real_t> y = test_data["y"].get<std::vector<real_t>>();
        std::vector<real_t> x_fit;
        std::vector<real_t> y_fit;
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
        auto y_eval_fit   = test_data["bspline_fit"]["y_eval"].get<std::vector<real_t>>();
        auto knots_fit    = test_data["bspline_fit"]["knots"].get<std::vector<real_t>>();
        auto ctrl_pts_fit = test_data["bspline_fit"]["ctrl"].get<std::vector<real_t>>();

        bspline.fit(x_fit, y_fit);

        REQUIRE_THAT(bspline.get_knots(), VectorsWithinAbsRel(knots_fit));

        if (BoundaryCondition::PERIODIC == BSplineType::boundary_condition_type)
        {
          SKIP("SciPy currently does not support fitting for periodic BSplines");
        }

        REQUIRE_THAT(bspline.get_control_points(), VectorsWithinAbsRel(ctrl_pts_fit));
        REQUIRE_THAT(bspline.evaluate(x_eval), VectorsWithinAbsRel(y_eval_fit));
      }

      SECTION("interpolate(...)")
      {

        std::vector<real_t> x_interp   = test_data["x"].get<std::vector<real_t>>();
        std::vector<real_t> y_interp   = test_data["y"].get<std::vector<real_t>>();
        auto [conds_left, conds_right] = test_data["conditions_interp"]
                                             .get<std::pair<
                                                 std::vector<std::pair<real_t, real_t>>,
                                                 std::vector<std::pair<real_t, real_t>>>>();
        std::vector<lsq::Condition<real_t>> additional_conditions{};
        additional_conditions.reserve(conds_left.size() + conds_right.size());
        for (auto [derivative_order, value] : conds_left)
        {
          additional_conditions.emplace_back(x_interp.front(), value, derivative_order);
        }
        for (auto [derivative_order, value] : conds_right)
        {
          additional_conditions.emplace_back(x_interp.back(), value, derivative_order);
        }

        auto y_eval_interp   = test_data["bspline_interp"]["y_eval"].get<std::vector<real_t>>();
        auto knots_interp    = test_data["bspline_interp"]["knots"].get<std::vector<real_t>>();
        auto ctrl_pts_interp = test_data["bspline_interp"]["ctrl"].get<std::vector<real_t>>();
        if constexpr (BoundaryCondition::OPEN == BSplineType::boundary_condition_type)
        {
          if constexpr (Curve::UNIFORM == BSplineType::curve_type)
          {
            SKIP("SciPy implementation forces clamped boundary condition, but we cannot add "
                 "explicit padding to uniform BSplines.");
          }
          x_interp.insert(x_interp.begin(), degree, x_interp.front());
          x_interp.insert(x_interp.end(), degree, x_interp.back());
          y_interp.insert(y_interp.begin(), degree, y_interp.front());
          y_interp.insert(y_interp.end(), degree, y_interp.back());
        }

        bspline.interpolate(x_interp, y_interp, additional_conditions);

        if (BoundaryCondition::PERIODIC == BSplineType::boundary_condition_type and
            (degree == 1 or degree % 2 == 0))
        {
          SKIP("SciPy implementation is similar only for odd degrees > 1");
        }

        REQUIRE_THAT(bspline.get_knots(), VectorsWithinAbsRel(knots_interp));
        REQUIRE_THAT(bspline.get_control_points(), VectorsWithinAbsRel(ctrl_pts_interp));
        REQUIRE_THAT(bspline.evaluate(x_eval), VectorsWithinAbsRel(y_eval_interp));
      }
    }
  }
}
