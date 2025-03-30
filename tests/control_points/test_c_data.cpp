// Standard includes

// Third-party includes
#include <catch2/catch_test_macros.hpp>
#include <vector>

// BSplineX includes
#include "BSplineX/control_points/c_data.hpp"

using namespace bsplinex;
using namespace bsplinex::control_points;

TEST_CASE("control_points::Data<double> data{data_vec}", "[c_data]")
{
  std::vector<double> const data_vec{0.0, 1.3, 2.2, 4.9, 13.2};
  Data<double> const data{data_vec};

  SECTION("data.size()") { REQUIRE(data.size() == 5); }
  SECTION("data.at(...)")
  {
    for (size_t i{0}; i < data.size(); i++)
    {
      REQUIRE(data.at(i) == data_vec.at(i));
    }
  }
}
