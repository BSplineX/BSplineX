// Third-party includes
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// BSplineX includes
#include "control_points/control_points.hpp"

using namespace Catch::Matchers;
using namespace bsplinex;
using namespace bsplinex::control_points;

TEST_CASE(
    "control_points::ControlPoints<T, BC> "
    "control_points{control_points::Data<T> data, degree}",
    "[control_points]"
)
{
  std::vector<double> data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  Data<double> data{data_vec};
  size_t degree{3};
  ControlPoints<double, BoundaryCondition::PERIODIC> control_points{data, degree};

  SECTION("control_points.size()") { REQUIRE(control_points.size() == data.size() + degree); }
  SECTION("control_points.at(...)")
  {
    for (size_t i{0}; i < data.size(); i++)
    {
      REQUIRE(control_points.at(i) == data.at(i));
    }
    REQUIRE(control_points.at(data.size()) == 0.1);
    REQUIRE(control_points.at(data.size() + 1) == 1.3);
    REQUIRE(control_points.at(data.size() + 2) == 2.2);
  }
}
