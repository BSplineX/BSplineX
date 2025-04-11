// Standard includes
#include <cstddef>
#include <vector>

// Third-party includes
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>

// BSplineX includes
#include "BSplineX/knots/knots.hpp"
#include "BSplineX/knots/t_data.hpp"
#include "BSplineX/types.hpp"
#include "matchers.hpp"

using namespace Catch::Matchers;
using namespace bsplinex;
using namespace bsplinex::knots;

TEST_CASE("knots::Knots<T, C, BC, EXT> knots{knots::Data<T, C> data, degree}", "[knots]")
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 2.2, 4.9, 6.3, 6.3, 6.3, 13.2};
  size_t const n{data_vec.size()};
  Data<double, Curve::NON_UNIFORM> const data{data_vec};
  size_t constexpr degree{3};
  Knots<double, Curve::NON_UNIFORM, BoundaryCondition::PERIODIC, Extrapolation::PERIODIC> const
      knots{data, degree};

  SECTION("knots.size()") { REQUIRE(knots.size() == data.size() + (2 * degree)); }
  SECTION("knots.at(...)")
  {
    REQUIRE_THAT(
        knots.at(0), WithinAbsRel(data_vec.at(0) - (data_vec.at(n - 1) - data_vec.at(n - 4)))
    );
    REQUIRE_THAT(
        knots.at(1), WithinAbsRel(data_vec.at(0) - (data_vec.at(n - 1) - data_vec.at(n - 3)))
    );
    REQUIRE_THAT(
        knots.at(2), WithinAbsRel(data_vec.at(0) - (data_vec.at(n - 1) - data_vec.at(n - 2)))
    );
    for (size_t i{0}; i < data.size(); i++)
    {
      REQUIRE(knots.at(i + degree) == data.at(i));
    }
    REQUIRE_THAT(
        knots.at(data.size() + degree),
        WithinAbsRel(data_vec.at(n - 1) + (data_vec.at(1) - data_vec.at(0)))
    );
    REQUIRE_THAT(
        knots.at(data.size() + degree + 1),
        WithinAbsRel(data_vec.at(n - 1) + (data_vec.at(2) - data_vec.at(0)))
    );
    REQUIRE_THAT(
        knots.at(data.size() + degree + 2),
        WithinAbsRel(data_vec.at(n - 1) + (data_vec.at(3) - data_vec.at(0)))
    );
  }
  // di
  SECTION("knots.find()")
  {
    std::vector<double> const values_find{-1.0, 0.1, 2.0, 2.2, 6.3, 13.2, 14.0};
    std::vector<size_t> const expected_indices{10, 3, 4, 6, 10, 10, 3};
    for (size_t i{0}; i < values_find.size(); i++)
    {
      REQUIRE(knots.find(values_find.at(i)).first == expected_indices.at(i));
    }
  }
}
