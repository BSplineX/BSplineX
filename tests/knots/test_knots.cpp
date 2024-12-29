// Third-party includes
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// BSplineX includes
#include "BSplineX/knots/knots.hpp"

using namespace Catch::Matchers;
using namespace bsplinex;
using namespace bsplinex::knots;

TEST_CASE("knots::Knots<T, C, BC, EXT> knots{knots::Data<T, C> data, degree}", "[knots]")
{
  std::vector<double> data_vec{0.1, 1.3, 2.2, 2.2, 4.9, 6.3, 6.3, 6.3, 13.2};
  Data<double, Curve::NON_UNIFORM> data{data_vec};
  size_t degree{3};
  Knots<double, Curve::NON_UNIFORM, BoundaryCondition::PERIODIC, Extrapolation::PERIODIC> knots{
      data, degree
  };

  SECTION("knots.size()") { REQUIRE(knots.size() == data.size() + 2 * degree); }
  SECTION("knots.at(...)")
  {
    REQUIRE_THAT(knots.at(0), WithinRel(-6.8));
    REQUIRE_THAT(knots.at(1), WithinRel(-6.8));
    REQUIRE_THAT(knots.at(2), WithinRel(-6.8));
    for (size_t i{0}; i < data.size(); i++)
    {
      REQUIRE(knots.at(i + degree) == data.at(i));
    }
    REQUIRE_THAT(knots.at(data.size() + degree), WithinRel(14.4));
    REQUIRE_THAT(knots.at(data.size() + degree + 1), WithinRel(15.3));
    REQUIRE_THAT(knots.at(data.size() + degree + 2), WithinRel(15.3));
  }
  SECTION("knots.find()")
  {
    REQUIRE(knots.find(-1.0).first == 10);
    REQUIRE(knots.find(0.1).first == 3);
    REQUIRE(knots.find(2.0).first == 4);
    REQUIRE(knots.find(2.2).first == 6);
    REQUIRE(knots.find(6.3).first == 10);
    REQUIRE(knots.find(13.2).first == 10);
    REQUIRE(knots.find(14.0).first == 3);
  }
}
