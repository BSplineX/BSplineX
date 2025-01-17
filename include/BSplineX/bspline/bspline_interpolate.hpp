#ifndef BSPLINE_INTERPOLATE_HPP
#define BSPLINE_INTERPOLATE_HPP

// Standard includes
#include <vector>

// BSplineX includes
#include "BSplineX/bspline/bspline_lsq.hpp"
#include "BSplineX/control_points/control_points.hpp"
#include "BSplineX/knots/knots.hpp"
#include "BSplineX/types.hpp"

namespace bsplinex::interpolate
{

template <typename T, Curve C, BoundaryCondition BC, Extrapolation EXT>
control_points::ControlPoints<T, BC> interpolate(
    size_t degree,
    knots::Knots<T, C, BC, EXT> const &knots,
    std::vector<T> const &y,
    std::vector<lsq::Condition<T>> const &additional_conditions
)
{
  assertm(
      additional_conditions.size() == degree - 1,
      "There must be exactly degree - 1 additional conditions."
  );

  size_t const &num_rows = knots.size() - 2 * degree;

  std::vector<T> x;
  x.reserve(num_rows);
  for (size_t i{0}; i < num_rows; i++)
  {
    x.push_back(knots.at(i + degree));
  }

  return lsq::lsq(degree, knots, x, y, additional_conditions);
}

} // namespace bsplinex::interpolate

#endif
