// Standard includes
#include <cstddef>
#include <stdexcept>
#include <vector>

// Third-party includes
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>

// BSplineX includes
#include "BSplineX/knots/t_data.hpp"
#include "BSplineX/knots/t_padder.hpp"
#include "BSplineX/types.hpp"
#include "matchers.hpp"

using namespace Catch::Matchers;
using namespace bsplinex;
using namespace bsplinex::knots;

TEST_CASE(
    "knots::Padder<T, C, BoundaryCondition::OPEN> padder{knots::Data<T, C> "
    "data, degree}",
    "[t_padder]"
)
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  Data<double, Curve::NON_UNIFORM> const data{data_vec};
  Padder<double, Curve::NON_UNIFORM, BoundaryCondition::OPEN> const padder{data, 3};

  SECTION("padder.size()") { REQUIRE(padder.size() == 0); }
  SECTION("padder.size_left()") { REQUIRE(padder.size_left() == 0); }
  SECTION("padder.size_right()") { REQUIRE(padder.size_right() == 0); }
  SECTION("padder.left()") { REQUIRE_THROWS_AS(padder.left(0), std::runtime_error); }
  SECTION("padder.right()") { REQUIRE_THROWS_AS(padder.right(0), std::runtime_error); }
}

TEST_CASE(
    "knots::Padder<T, C, BoundaryCondition::CLAMPED> padder{knots::Data<T, C> "
    "data, degree}",
    "[t_padder]"
)
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  Data<double, Curve::NON_UNIFORM> const data{data_vec};
  Padder<double, Curve::NON_UNIFORM, BoundaryCondition::CLAMPED> const padder{data, 3};

  SECTION("padder.size()") { REQUIRE(padder.size() == 6); }
  SECTION("padder.size_left()") { REQUIRE(padder.size_left() == 3); }
  SECTION("padder.size_right()") { REQUIRE(padder.size_right() == 3); }
  SECTION("padder.left(...)")
  {
    for (size_t i{0}; i < padder.size_left(); i++)
    {
      REQUIRE(padder.left(i) == data.at(0));
    }
  }
  SECTION("padder.right(...)")
  {
    for (size_t i{0}; i < padder.size_right(); i++)
    {
      REQUIRE(padder.right(i) == data.at(data.size() - 1));
    }
  }
}

TEST_CASE(
    "knots::Padder<T, C, BoundaryCondition::PERIODIC> padder{knots::Data<T, C> "
    "data, degree}",
    "[t_padder]"
)
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  size_t const n = data_vec.size();
  Data<double, Curve::NON_UNIFORM> const data{data_vec};
  Padder<double, Curve::NON_UNIFORM, BoundaryCondition::PERIODIC> const padder{data, 3};

  SECTION("padder.size()") { REQUIRE(padder.size() == 6); }
  SECTION("padder.size_left()") { REQUIRE(padder.size_left() == 3); }
  SECTION("padder.size_right()") { REQUIRE(padder.size_right() == 3); }
  SECTION("padder.left(...)")
  {
    REQUIRE_THAT(
        padder.left(0), WithinAbsRel(data_vec.at(0) - (data_vec.at(n - 1) - data_vec.at(n - 4)))
    );
    REQUIRE_THAT(
        padder.left(1), WithinAbsRel(data_vec.at(0) - (data_vec.at(n - 1) - data_vec.at(n - 3)))
    );
    REQUIRE_THAT(
        padder.left(2), WithinAbsRel(data_vec.at(0) - (data_vec.at(n - 1) - data_vec.at(n - 2)))
    );
  }
  SECTION("padder.right(...)")
  {
    REQUIRE_THAT(
        padder.right(0), WithinAbsRel(data_vec.at(n - 1) + (data_vec.at(1) - data_vec.at(0)))
    );
    REQUIRE_THAT(
        padder.right(1), WithinAbsRel(data_vec.at(n - 1) + (data_vec.at(2) - data_vec.at(0)))
    );
    REQUIRE_THAT(
        padder.right(2), WithinAbsRel(data_vec.at(n - 1) + (data_vec.at(3) - data_vec.at(0)))
    );
  }
}
