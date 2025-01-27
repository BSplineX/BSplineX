#ifndef BSPLINEX_BSPLINE_BSPLINE_FACTORY_PERIODIC_HPP
#define BSPLINEX_BSPLINE_BSPLINE_FACTORY_PERIODIC_HPP

// Standard
#include <vector>

// BSplineX
#include "BSplineX/bspline/bspline_types.hpp"

namespace bsplinex::factory
{

template <typename T = double>
inline types::PeriodicUniform<T>
periodic_uniform(size_t degree, T begin, T end, size_t num_elems, std::vector<T> const &ctrl_points)
{
  return types::PeriodicUniform<T>{{begin, end, num_elems}, {ctrl_points}, degree};
}

template <typename T = double>
inline types::PeriodicUniform<T> periodic_uniform(size_t degree, T begin, T end, size_t num_elems)
{
  return periodic_uniform<T>(degree, begin, end, num_elems, std::vector<T>(num_elems - 1));
}

template <typename T = double>
inline types::PeriodicUniform<T> periodic_uniform(size_t degree)
{
  constexpr T begin{(T)0};
  constexpr T end{(T)1};
  size_t const num_elems{2 * degree};
  return periodic_uniform<T>(degree, begin, end, num_elems, std::vector<T>(num_elems - 1));
}

template <typename T = double>
inline types::PeriodicNonUniform<T>
periodic_nonuniform(size_t degree, std::vector<T> const &knots, std::vector<T> const &ctrl_points)
{
  return types::PeriodicNonUniform<T>{{knots}, {ctrl_points}, degree};
}

template <typename T = double>
inline types::PeriodicNonUniform<T> periodic_nonuniform(size_t degree, std::vector<T> const &knots)
{
  return periodic_nonuniform<T>(degree, knots, std::vector<T>(knots.size() - 1));
}

template <typename T = double>
inline types::PeriodicNonUniform<T> periodic_nonuniform(size_t degree)
{
  std::vector<T> const knots(2 * degree);
  return periodic_nonuniform<T>(degree, knots, std::vector<T>(knots.size() - 1));
}

} // namespace bsplinex::factory

#endif
