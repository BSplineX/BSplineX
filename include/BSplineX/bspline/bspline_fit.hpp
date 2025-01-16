#ifndef BSPLINE_FIT_HPP
#define BSPLINE_FIT_HPP

// Standard includes
#include <vector>

// BSplineX includes
#include "BSplineX/bspline/bspline_lsq.hpp"
#include "BSplineX/control_points/control_points.hpp"
#include "BSplineX/knots/knots.hpp"
#include "BSplineX/types.hpp"

namespace bsplinex::fit
{

template <typename T, Curve C, BoundaryCondition BC, Extrapolation EXT>
control_points::ControlPoints<T, BC>
fit(size_t degree,
    knots::Knots<T, C, BC, EXT> const &knots,
    std::vector<T> const &x,
    std::vector<T> const &y)
{
  return lsq::lsq(degree, knots, x, y, {});
}

} // namespace bsplinex::fit

#endif
