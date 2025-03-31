// Standard includes
#include <stdexcept>
#include <vector>

// Third-party includes
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// BSplineX includes
#include "BSplineX/control_points/c_data.hpp"
#include "BSplineX/control_points/c_padder.hpp"
#include "BSplineX/types.hpp"

using namespace Catch::Matchers;
using namespace bsplinex;
using namespace bsplinex::control_points;

TEST_CASE(
    "control_points::Padder<T, BC> padder{control_points::Data<T> data, degree}", "[c_padder]"
)
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  Data<double> data{data_vec};
  Padder<double, BoundaryCondition::OPEN> const padder{data, 3};

  SECTION("padder.size()") { REQUIRE(padder.size() == 0); }
  SECTION("padder.size_right()") { REQUIRE(padder.size_right() == 0); }
  SECTION("padder.right()") { REQUIRE_THROWS_AS(padder.right(0), std::runtime_error); }
}

TEST_CASE(
    "control_points::Padder<T, BoundaryCondition::PERIODIC>"
    "padder{control_points::Data<T> data, degree}",
    "[c_padder]"
)
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  Data<double> data{data_vec};
  Padder<double, BoundaryCondition::PERIODIC> const padder{data, 3};

  SECTION("padder.size()") { REQUIRE(padder.size() == 3); }
  SECTION("padder.size_right()") { REQUIRE(padder.size_right() == 3); }
  SECTION("padder.right(...)")
  {
    REQUIRE_THAT(padder.right(0), WithinAbs(data_vec[0], 0));
    REQUIRE_THAT(padder.right(1), WithinAbs(data_vec[1], 0));
    REQUIRE_THAT(padder.right(2), WithinAbs(data_vec[2], 0));
  }
}
