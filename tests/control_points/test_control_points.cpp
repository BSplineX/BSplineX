// Standard includes
#include <cstddef>
#include <vector>

// Third-party includes
#include <catch2/catch_test_macros.hpp>

// BSplineX includes
#include "BSplineX/control_points/c_data.hpp"
#include "BSplineX/control_points/control_points.hpp"
#include "BSplineX/types.hpp"

using namespace bsplinex;
using namespace bsplinex::control_points;

TEST_CASE(
    "control_points::ControlPoints<T, BC> "
    "control_points{control_points::Data<T> data, degree}",
    "[control_points]"
)
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  Data<double> const data{data_vec};
  size_t constexpr degree{3};
  ControlPoints<double, BoundaryCondition::PERIODIC> const control_points{data, degree};

  SECTION("control_points.size()") { REQUIRE(control_points.size() == data.size() + degree); }
  SECTION("control_points.at(...)")
  {
    for (size_t i{0}; i < data.size(); i++)
    {
      REQUIRE(control_points.at(i) == data.at(i));
    }
    REQUIRE(control_points.at(data.size()) == data.at(0));
    REQUIRE(control_points.at(data.size() + 1) == data.at(1));
    REQUIRE(control_points.at(data.size() + 2) == data.at(2));
  }
}
