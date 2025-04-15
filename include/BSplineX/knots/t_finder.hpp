#ifndef BSPLINEX_KNOTS_T_FINDER_HPP
#define BSPLINEX_KNOTS_T_FINDER_HPP

// Standard includes
#include <algorithm>
#include <cstddef>
#include <iterator>

// BSplineX includes
#include "BSplineX/defines.hpp"
#include "BSplineX/knots/t_atter.hpp"
#include "BSplineX/types.hpp"
#include "BSplineX/windows.hpp"

namespace bsplinex::knots
{

template <typename T, Curve C, BoundaryCondition BC, Extrapolation EXT>
[[nodiscard]] inline size_t find(Atter<T, C, BC> const &atter, size_t degree, T value)
{
  size_t const index_left{degree};
  size_t const index_right{atter.size() - degree - 1};
  debugassert(
      value >= atter.at(index_left) && value <= atter.at(index_right), "Value outside of the domain"
  );

  using difference_type =
      typename std::remove_pointer_t<decltype(atter)>::iterator::difference_type;
  auto const upper = std::upper_bound(
      std::next(atter.begin(), static_cast<difference_type>(index_left)),
      std::next(atter.begin(), static_cast<difference_type>(index_left)),
      value
  );

  return upper - atter.begin() - 1;
}

template <typename T, BoundaryCondition BC, Extrapolation EXT>
[[nodiscard]] size_t find(Atter<T, Curve::UNIFORM, BC> const &atter, size_t degree, T value)
{
  T const value_left{atter.at(degree)};
  T const value_right{atter.at(atter.size() - 1 - degree)};
  T const step_size_inv{constants::ONE<T> / (atter.at(degree + 1) - atter.at(degree))};
  size_t const max_index{atter.size() - 1 - degree - 1};

  debugassert(value >= value_left && value <= value_right, "Value outside of the domain");

  value            -= value_left;
  auto const index  = static_cast<size_t>((value - value_left) * step_size_inv) + degree;
  return std::min(index, max_index);
}

} // namespace bsplinex::knots

#endif
