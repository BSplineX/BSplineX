// Standard includes
#include <cstddef>
#include <stdexcept>
#include <vector>

// Third-party includes
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>

// BSplineX includes
#include "BSplineX/knots/t_atter.hpp"
#include "BSplineX/knots/t_data.hpp"
#include "BSplineX/knots/t_extrapolator.hpp"
#include "BSplineX/types.hpp"
#include "matchers.hpp"

using namespace Catch::Matchers;
using namespace bsplinex;
using namespace bsplinex::knots;

TEST_CASE(
    "knots::Extrapolator<T, C, BC, Extrapolation::NONE> extrapolator{atter}", "[t_extrapolator]"
)
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  Data<double, Curve::NON_UNIFORM> const data{data_vec};
  size_t constexpr degree{3};
  Atter<double, Curve::NON_UNIFORM, BoundaryCondition::CLAMPED> const atter{data, degree};
  Extrapolator<double, Curve::NON_UNIFORM, BoundaryCondition::CLAMPED, Extrapolation::NONE> const
      extrapolator{atter, degree};

  SECTION("extrapolator.extrapolate()")
  {
    double constexpr out_left{-1.0};
    double constexpr out_right{15.0};
    REQUIRE_THROWS_AS(extrapolator.extrapolate(out_left), std::runtime_error);
    REQUIRE_THROWS_AS(extrapolator.extrapolate(out_right), std::runtime_error);
  }
}

TEST_CASE(
    "knots::Extrapolator<T, C, BC, Extrapolation::CONSTANT> "
    "extrapolator{atter}",
    "[t_extrapolator]"
)
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  Data<double, Curve::NON_UNIFORM> const data{data_vec};
  size_t constexpr degree{3};
  Atter<double, Curve::NON_UNIFORM, BoundaryCondition::CLAMPED> const atter{data, degree};
  Extrapolator<
      double,
      Curve::NON_UNIFORM,
      BoundaryCondition::CLAMPED,
      Extrapolation::CONSTANT> const extrapolator{atter, degree};

  SECTION("extrapolator.extrapolate()")
  {
    double constexpr out_left{-1.0};
    double constexpr out_right{15.0};
    REQUIRE(extrapolator.extrapolate(out_left) == atter.at(0));
    REQUIRE(extrapolator.extrapolate(out_right) == atter.at(atter.size() - 1));
  }
}

TEST_CASE(
    "knots::Extrapolator<T, C, BC, Extrapolation::PERIODIC> "
    "extrapolator{atter}",
    "[t_extrapolator]"
)
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  Data<double, Curve::NON_UNIFORM> const data{data_vec};
  size_t constexpr degree{3};
  Atter<double, Curve::NON_UNIFORM, BoundaryCondition::PERIODIC> const atter{data, degree};
  Extrapolator<
      double,
      Curve::NON_UNIFORM,
      BoundaryCondition::PERIODIC,
      Extrapolation::PERIODIC> const extrapolator{atter, degree};
  double const period = data.at(data.size() - 1) - data.at(0);

  SECTION("extrapolator.extrapolate()")
  {
    double constexpr out_left{-1.0};
    double constexpr out_right{15.0};
    REQUIRE_THAT(
        extrapolator.extrapolate(out_left - (2 * period)), WithinAbsRel(out_left + period)
    );
    REQUIRE_THAT(extrapolator.extrapolate(out_left - period), WithinAbsRel(out_left + period));
    REQUIRE_THAT(extrapolator.extrapolate(out_left), WithinAbsRel(out_left + period));
    REQUIRE_THAT(extrapolator.extrapolate(out_right), WithinAbsRel(out_right - period));
    REQUIRE_THAT(extrapolator.extrapolate(out_right + period), WithinAbsRel(out_right - period));
    REQUIRE_THAT(
        extrapolator.extrapolate(out_right + (2 * period)), WithinAbsRel(out_right - period)
    );
  }
}
