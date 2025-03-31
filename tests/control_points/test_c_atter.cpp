// Standard includes
#include <cstddef>
#include <vector>

// Third-party includes
#include <catch2/catch_test_macros.hpp>

// BSplineX includes
#include "BSplineX/control_points/c_atter.hpp"
#include "BSplineX/control_points/c_data.hpp"
#include "BSplineX/types.hpp"

using namespace bsplinex;
using namespace bsplinex::control_points;

TEST_CASE("control_points::Atter<T, BC> atter{control_points::Data<T> data, degree}", "[c_atter]")
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  Data<double> const data{data_vec};
  size_t constexpr degree{3};
  Atter<double, BoundaryCondition::PERIODIC> const atter{data, degree};

  SECTION("atter.size()") { REQUIRE(atter.size() == data.size() + degree); }
  SECTION("atter.at(...)")
  {
    for (size_t i{0}; i < data.size(); i++)
    {
      REQUIRE(atter.at(i) == data.at(i));
    }
    REQUIRE(atter.at(data.size()) == data.at(0));
    REQUIRE(atter.at(data.size() + 1) == data.at(1));
    REQUIRE(atter.at(data.size() + 2) == data.at(2));
  }
}
