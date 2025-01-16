#ifndef BSPLINE_INTERPOLATE_HPP
#define BSPLINE_INTERPOLATE_HPP

// Standard includes
#include <algorithm>
#include <functional>
#include <utility>
#include <vector>

// BSplineX includes
#include "BSplineX/bspline/bspline_lsq.hpp"
#include "BSplineX/control_points/control_points.hpp"
#include "BSplineX/knots/knots.hpp"
#include "BSplineX/types.hpp"

namespace bsplinex::interpolate
{

template <typename T, BoundaryCondition BC, Extrapolation EXT>
std::pair<knots::Knots<T, Curve::UNIFORM, BC, EXT>, control_points::ControlPoints<T, BC>> uniform(
    size_t degree,
    T begin,
    T end,
    size_t num_elems,
    std::vector<T> const &y,
    std::vector<lsq::Condition<T>> const &additional_conditions
)
{
  size_t step_size = (end - begin) / (num_elems - 1);
  knots::Knots<T, Curve::UNIFORM, BC, EXT> knots{};
  if constexpr (BoundaryCondition::OPEN == BC)
  {
    new (&knots) knots::Knots<T, Curve::UNIFORM, BC, EXT>{
        {begin - degree * step_size, end + degree * step_size, num_elems + 2 * degree}, degree
    };
  }
  else
  {
    new (&knots) knots::Knots<T, Curve::UNIFORM, BC, EXT>{{begin, end, num_elems}, degree};
  }

  std::vector<T> x;
  x.reserve(y.size());
  for (size_t i{0}; i < x.size(); i++)
  {
    x.push_back(begin + i * step_size);
  }

  control_points::ControlPoints<T, BC> ctrl_pts =
      lsq::lsq(degree, knots, x, y, additional_conditions);

  return std::make_pair(std::move(knots), std::move(ctrl_pts));
}

template <typename T, BoundaryCondition BC, Extrapolation EXT>
std::pair<knots::Knots<T, Curve::NON_UNIFORM, BC, EXT>, control_points::ControlPoints<T, BC>>
non_uniform(
    size_t degree,
    std::vector<T> const &x,
    std::vector<T> const &y,
    [[maybe_unused]] std::vector<T> const &padding,
    std::vector<lsq::Condition<T>> const &additional_conditions
)
{
  // NOTE: thank the STL for this wonderful backwards built sort check. Think it as if std::less
  // is <= and std::less_equal is <.
  if (!std::is_sorted(x.begin(), x.end(), std::less_equal<T>{}))
  {
    throw std::runtime_error("x must be sorted w.r.t. operator <");
  }

  knots::Knots<T, Curve::NON_UNIFORM, BC, EXT> knots{};
  if constexpr (BoundaryCondition::OPEN == BC)
  {
    if (padding.size() != 2 * degree)
    {
      throw std::runtime_error("padding must be = 2 * degree");
    }
    size_t open_knots_size = x.size() + 2 * degree;
    std::vector<T> open_knots{};
    open_knots.reserve(open_knots_size);
    for (size_t i{0}; i < degree; i++)
    {
      open_knots.push_back(padding.at(i));
    }
    for (auto const &elem : x)
    {
      open_knots.push_back(elem);
    }
    for (size_t i{0}; i < degree; i++)
    {
      open_knots.push_back(padding.at(i + degree));
    }

    new (&knots) knots::Knots<T, Curve::NON_UNIFORM, BC, EXT>{{open_knots}, degree};
  }
  else
  {
    new (&knots) knots::Knots<T, Curve::NON_UNIFORM, BC, EXT>{{x}, degree};
  }

  control_points::ControlPoints<T, BC> ctrl_pts =
      lsq::lsq(degree, knots, x, y, additional_conditions);

  return std::make_pair(std::move(knots), std::move(ctrl_pts));
}

} // namespace bsplinex::interpolate

#endif
