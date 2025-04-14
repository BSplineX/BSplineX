// Standard includes
#include <cstddef>
#include <vector>

// Third-party includes
#include <catch2/catch_test_macros.hpp>

// BSplineX includes
#include "BSplineX/knots/t_atter.hpp"
#include "BSplineX/knots/t_data.hpp"
#include "BSplineX/knots/t_finder.hpp"
#include "BSplineX/types.hpp"

using namespace bsplinex;
using namespace bsplinex::knots;

TEST_CASE(
    "knots::Finder<double, NON_UNIFORM, PERIODIC, PERIODIC> "
    "finder{atter}",
    "[t_finder]"
)
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 2.2, 4.9, 6.3, 6.3, 6.3, 13.2};
  Data<double, Curve::NON_UNIFORM> const data{data_vec};
  size_t constexpr degree{3};
  Atter<double, Curve::NON_UNIFORM, BoundaryCondition::PERIODIC> const atter{data, degree};

  SECTION("finder.find()")
  {
    std::vector<double> const values_find{0.1, 2.0, 2.2, 6.3};
    std::vector<size_t> const expected_indices{degree, degree + 1, degree + 3, degree + 7};
    for (size_t i{0}; i < values_find.size(); i++)
    {
      REQUIRE(
          knots::find<
              double,
              Curve::NON_UNIFORM,
              BoundaryCondition::PERIODIC,
              Extrapolation::PERIODIC>(atter, degree, values_find.at(i)) == expected_indices.at(i)
      );
    }
  }
}
