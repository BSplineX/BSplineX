// Standard includes
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <tuple>
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
#include "BSplineX/control_points/c_data.hpp"
#include "BSplineX/knots/t_data.hpp"
#include "BSplineX/types.hpp"
#include "BSplineX/windows.hpp"
#include "matchers.hpp"

using real_t = double;

using namespace Catch::Matchers;
using namespace bsplinex;
using namespace bsplinex::bspline;

// clang-format off
#define BSPLINE_TEST_TYPES \
types::OpenUniform<real_t>, \
types::OpenUniformConstant<real_t>, \
types::OpenNonUniform<real_t>, \
types::OpenNonUniformConstant<real_t>, \
types::ClampedUniform<real_t>, \
types::ClampedUniformConstant<real_t>, \
types::ClampedNonUniform<real_t>, \
types::ClampedNonUniformConstant<real_t>, \
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
  return path.string();
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

    if (BoundaryCondition::CLAMPED == BSplineType::boundary_condition_type or
        BoundaryCondition::PERIODIC == BSplineType::boundary_condition_type)
    {
      knots_begin = knots.at(degree);
      knots_end   = knots.at(knots.size() - degree - 1);
      num_knots   = knots.size() - 2 * degree;
    }

    return BSplineType(
        knots::Data<real_t, BSplineType::curve_type>{knots_begin, knots_end, num_knots},
        control_points::Data<real_t>{ctrl_pts},
        degree
    );
  }
  else
  {
    static_assert(Curve::NON_UNIFORM == BSplineType::curve_type);
    if (BoundaryCondition::CLAMPED == BSplineType::boundary_condition_type or
        BoundaryCondition::PERIODIC == BSplineType::boundary_condition_type)
    {
      using difference_type = typename std::vector<real_t>::iterator::difference_type;

      auto const deg = static_cast<difference_type>(degree);
      knots = std::vector<real_t>(std::next(knots.begin(), deg), std::prev(knots.end(), deg));
    }
    return BSplineType(
        knots::Data<real_t, BSplineType::curve_type>{knots},
        control_points::Data<real_t>{ctrl_pts},
        degree
    );
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

      auto const domain = test_data["bspline"]["domain"].get<std::pair<real_t, real_t>>();
      SECTION("domain()")
      {
        auto [domain_left, domain_right] = bspline.domain();
        REQUIRE(domain_left == domain.first);
        REQUIRE(domain_right == domain.second);
      }
      SECTION("evaluate(...)")
      {
        auto derivative_data = test_data["bspline"];
        for (size_t derivative_order{0}; derivative_order <= degree; derivative_order++)
        {
          SECTION("evaluate(..., derivative_order=" + std::to_string(derivative_order) + ")")
          {
            auto const y_eval = derivative_data["y_eval"].get<std::vector<real_t>>();
            for (size_t i{0}; i < x_eval.size(); i++)
            {
              REQUIRE_THAT(
                  bspline.evaluate(x_eval.at(i), derivative_order), WithinAbsRel(y_eval.at(i))
              );
            }
            REQUIRE_THAT(bspline.evaluate(x_eval, derivative_order), VectorsWithinAbsRel(y_eval));
            SECTION("extrapolation")
            {
              if constexpr (Extrapolation::NONE == BSplineType::extrapolation_type)
              {
                REQUIRE_THROWS_WITH(
                    bspline.evaluate(domain.first - 1, derivative_order),
                    "Extrapolation explicitly set to NONE"
                );
                REQUIRE_THROWS_WITH(
                    bspline.evaluate(domain.second + 1, derivative_order),
                    "Extrapolation explicitly set to NONE"
                );
              }
              else if constexpr (Extrapolation::CONSTANT == BSplineType::extrapolation_type)
              {
                real_t const y_left  = bspline.evaluate(domain.first, derivative_order);
                real_t const y_right = bspline.evaluate(domain.second, derivative_order);
                constexpr size_t extrapolation_elems{50};
                for (size_t i{0}; i < extrapolation_elems; i++)
                {
                  REQUIRE_THAT(
                      bspline.evaluate(domain.first - static_cast<real_t>(i), derivative_order),
                      WithinAbsRel(y_left)
                  );
                  REQUIRE_THAT(
                      bspline.evaluate(domain.second + static_cast<real_t>(i), derivative_order),
                      WithinAbsRel(y_right)
                  );
                }
              }
              else if constexpr (Extrapolation::PERIODIC == BSplineType::extrapolation_type)
              {
                if (derivative_order < degree)
                {
                  REQUIRE_THAT(
                      bspline.evaluate(domain.first, derivative_order),
                      WithinAbsRel(bspline.evaluate(domain.second, derivative_order))
                  );
                }
                real_t const period = domain.second - domain.first;
                constexpr real_t periodic_extrapolation_tol{
                    1e-8
                }; // HACK: this may be correct, but requires an in-depth analysis
                for (real_t const x : x_eval)
                {
                  real_t const y = bspline.evaluate(x, derivative_order);
                  for (size_t p{1}; p <= 3; ++p)
                  {
                    real_t const delta   = static_cast<real_t>(p) * period;
                    real_t const x_right = x + delta;
                    real_t const x_left  = x - delta;
                    REQUIRE_THAT(
                        bspline.evaluate(x_right, derivative_order),
                        WithinAbsRel(y, periodic_extrapolation_tol, periodic_extrapolation_tol)
                    );
                    REQUIRE_THAT(
                        bspline.evaluate(x_left, derivative_order),
                        WithinAbsRel(y, periodic_extrapolation_tol, periodic_extrapolation_tol)
                    );
                  }
                }
              }
            }
          }
          derivative_data = derivative_data["derivative"];
        }
        SECTION("evaluate(..., invalid derivative_order)")
        {
          REQUIRE_THROWS_WITH(
              bspline.evaluate(x_eval.at(0), degree + 1), "derivative_order must be in [1, degree]"
          );
          REQUIRE_THROWS_WITH(
              bspline.evaluate(x_eval, degree + 1), "derivative_order must be in [1, degree]"
          );
        }
      }

      SECTION("derivative(...)")
      {
        auto derivative_data = test_data["bspline"]["derivative"];
        for (size_t derivative_order{1}; derivative_order <= degree; derivative_order++)
        {
          SECTION("derivative(derivative_order=" + std::to_string(derivative_order) + ")")
          {
            BSplineType const derivative = bspline.derivative(derivative_order);
            auto const y_eval_d          = derivative_data["y_eval"].get<std::vector<real_t>>();
            REQUIRE(derivative.get_degree() == derivative_data["degree"].get<size_t>());
            REQUIRE_THAT(derivative.evaluate(x_eval), VectorsWithinAbsRel(y_eval_d));
          }
          derivative_data = derivative_data["derivative"];
        }
        SECTION("derivative(invalid derivative_order)")
        {
          REQUIRE_THROWS_WITH(bspline.derivative(0), "derivative_order must be in [1, degree]");
          REQUIRE_THROWS_WITH(
              bspline.derivative(degree + 1), "derivative_order must be in [1, degree]"
          );
        }
      }

      SECTION("nnz_basis(...)")
      {
        for (size_t derivative_order = 0; derivative_order <= degree; ++derivative_order)
        {
          SECTION("nnz_basis(..., derivative_order=" + std::to_string(derivative_order) + ")")
          {
            std::vector<real_t> x_nnz;
            x_nnz.reserve(x_eval.size());
            std::copy_if(
                x_eval.begin(),
                x_eval.end(),
                std::back_inserter(x_nnz),
                [domain](real_t x) { return domain.first <= x and x <= domain.second; }
            );
            auto const &nnz_basis_data = test_data["bspline"]["nnz_basis"][derivative_order];
            REQUIRE(x_nnz.size() == nnz_basis_data.size());

            for (size_t i = 0; i < x_nnz.size(); ++i)
            {
              auto const [ref_index, ref_nnz_basis] =
                  nnz_basis_data[i].get<std::pair<size_t, std::vector<real_t>>>();
              auto const [index, nnz_basis] = bspline.nnz_basis(x_nnz.at(i), derivative_order);
              REQUIRE(index == ref_index);
              REQUIRE_THAT(nnz_basis, VectorsWithinAbsRel(ref_nnz_basis));
            }
          }
        }
      }

      SECTION("fit(...)")
      {
        auto const x = test_data["x"].get<std::vector<real_t>>();
        auto const y = test_data["y"].get<std::vector<real_t>>();
        std::vector<real_t> x_fit;
        std::vector<real_t> y_fit;
        x_fit.reserve(x.size());
        y_fit.reserve(y.size());
        for (size_t i{0}; i < x.size(); i++)
        {
          if (domain.first <= x.at(i) and x.at(i) <= domain.second)
          {
            x_fit.push_back(x.at(i));
            y_fit.push_back(y.at(i));
          }
        }

        bspline.fit(x_fit, y_fit);
        // test fitting
        constexpr real_t fit_points_tol{1e-2};

        // ----------------- basic tests -----------------
        REQUIRE_THAT(
            bspline.evaluate(x_fit), VectorsWithinAbsRel(y_fit, fit_points_tol, fit_points_tol)
        );
        // Periodic BSpline and k-1 derivatives must match at the extremes
        for (size_t derivative_order{0}; derivative_order < degree; derivative_order++)
        {
          if constexpr (BoundaryCondition::PERIODIC == BSplineType::boundary_condition_type)
          {
            REQUIRE_THAT(
                bspline.evaluate(domain.first, derivative_order),
                WithinAbsRel(bspline.evaluate(domain.second, derivative_order))
            );
          }
        }

        // ----------------- test against reference data -----------------
        auto const &bspline_data = test_data["bspline_fit"];

        auto const knots_fit = bspline_data["knots"].get<std::vector<real_t>>();
        REQUIRE_THAT(bspline.get_knots(), VectorsWithinAbsRel(knots_fit));
        if (BoundaryCondition::PERIODIC == BSplineType::boundary_condition_type)
        {
          SKIP("SciPy currently does not support fitting for periodic BSplines");
        }

        auto const ctrl_pts_fit = bspline_data["ctrl"].get<std::vector<real_t>>();
        REQUIRE_THAT(bspline.get_control_points(), VectorsWithinAbsRel(ctrl_pts_fit));

        auto derivative_data = bspline_data;
        for (size_t derivative_order = 0; derivative_order <= degree; ++derivative_order)
        {
          SECTION("fit+evaluate(...,  derivative_order=" + std::to_string(derivative_order) + ")")
          {
            auto const y_eval_fit     = derivative_data["y_eval"].get<std::vector<real_t>>();
            real_t const fit_data_tol = derivative_order == 0 ? RTOL<real_t> : 1e-5;
            REQUIRE_THAT(
                bspline.evaluate(x_eval, derivative_order),
                VectorsWithinAbsRel(y_eval_fit, fit_data_tol, fit_data_tol)
            );
          }
          derivative_data = derivative_data["derivative"];
        }
      }

      SECTION("interpolate(...)")
      {

        auto x_interp                        = test_data["x"].get<std::vector<real_t>>();
        auto y_interp                        = test_data["y"].get<std::vector<real_t>>();
        auto const [conds_left, conds_right] = test_data["conditions_interp"]
                                                   .get<std::pair<
                                                       std::vector<std::pair<size_t, real_t>>,
                                                       std::vector<std::pair<size_t, real_t>>>>();
        std::vector<lsq::InterpolantCondition<real_t>> additional_conditions{};
        additional_conditions.reserve(conds_left.size() + conds_right.size());
        for (auto const &[derivative_order, value] : conds_left)
        {
          additional_conditions.emplace_back(x_interp.front(), value, derivative_order);
        }
        for (auto const &[derivative_order, value] : conds_right)
        {
          additional_conditions.emplace_back(x_interp.back(), value, derivative_order);
        }

        if constexpr (BoundaryCondition::OPEN == BSplineType::boundary_condition_type)
        {
          if constexpr (Curve::UNIFORM == BSplineType::curve_type)
          {
            SKIP(
                "SciPy implementation forces clamped boundary condition, but we cannot add "
                "explicit padding to uniform BSplines."
            );
          }
          x_interp.insert(x_interp.begin(), degree, x_interp.front());
          x_interp.insert(x_interp.end(), degree, x_interp.back());
          y_interp.insert(y_interp.begin(), degree, y_interp.front());
          y_interp.insert(y_interp.end(), degree, y_interp.back());
        }

        bspline.interpolate(x_interp, y_interp, additional_conditions);

        // --------------- basic tests -----------------
        REQUIRE_THAT(bspline.evaluate(x_interp), VectorsWithinAbsRel(y_interp));

        if (BoundaryCondition::PERIODIC == BSplineType::boundary_condition_type and
            (degree == 1 or degree % 2 == 0))
        {
          SKIP("SciPy implementation is similar only for odd degrees > 1");
        }

        // --------------- test against reference data -----------------
        auto derivative_data = test_data["bspline_interp"];

        for (size_t derivative_order = 0; derivative_order <= degree; ++derivative_order)
        {
          SECTION(
              "interpolate+evaluate(..., derivative_order=" + std::to_string(derivative_order) + ")"
          )
          {
            BSplineType const derivative =
                derivative_order == 0 ? bspline : bspline.derivative(derivative_order);

            auto const knots_interp =
                test_data["bspline_interp"]["knots"].get<std::vector<real_t>>();
            auto const ctrl_pts_interp =
                test_data["bspline_interp"]["ctrl"].get<std::vector<real_t>>();
            auto const y_eval_interp = derivative_data["y_eval"].get<std::vector<real_t>>();

            REQUIRE_THAT(bspline.get_knots(), VectorsWithinAbsRel(knots_interp));
            REQUIRE_THAT(bspline.get_control_points(), VectorsWithinAbsRel(ctrl_pts_interp));
            real_t const interp_data_tol = derivative_order == 0 ? RTOL<real_t> : 1e-5;
            REQUIRE_THAT(
                derivative.evaluate(x_eval),
                VectorsWithinAbsRel(y_eval_interp, interp_data_tol, interp_data_tol)
            );

            // Periodic BSpline and k-1 derivatives must match at the extremes
            if (BoundaryCondition::PERIODIC == BSplineType::boundary_condition_type and
                derivative_order < degree)
            {
              REQUIRE_THAT(
                  derivative.evaluate(domain.first),
                  WithinAbsRel(derivative.evaluate(domain.second))
              );
            }
          }
          derivative_data = derivative_data["derivative"];
        }
      }

      SECTION("Constructors and assignments")
      {
        auto const bspline_base = build_bspline<BSplineType>(
            {test_data["bspline"]["knots"].get<std::vector<real_t>>()},
            {test_data["bspline"]["ctrl"].get<std::vector<real_t>>()},
            degree
        );
        std::ignore =
            bspline_base.evaluate(domain.first, degree); // force computation of derivatives
        auto require_equals = [&bspline_base](BSplineType const &actual)
        {
          REQUIRE(actual == bspline_base);
          for (size_t derivative_order{1}; derivative_order <= bspline_base.get_degree();
               derivative_order++)
          {
            REQUIRE(
                actual.derivative(derivative_order) == bspline_base.derivative(derivative_order)
            );
          }
        };
        SECTION("Copy constructor")
        {
          BSplineType const new_bspline{
              bspline
          }; // NOLINT(performance-unnecessary-copy-initialization)
          require_equals(new_bspline);
        }

        SECTION("Copy assignment operator")
        {
          BSplineType new_bspline{};
          new_bspline = bspline;
          require_equals(new_bspline);
        }

        SECTION("Move constructor")
        {
          BSplineType tmp_bspline{bspline};
          BSplineType const new_bspline{std::move(tmp_bspline)};
          require_equals(new_bspline);
        }

        SECTION("Move assignment operator")
        {
          BSplineType tmp_bspline{bspline};
          BSplineType new_bspline{};
          new_bspline = std::move(tmp_bspline);
          require_equals(new_bspline);
        }
      }
    }
  }
}
