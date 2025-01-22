#ifndef BSPLINE_FACTORY_OPEN_HPP
#define BSPLINE_FACTORY_OPEN_HPP

// Standard
#include <vector>

// BSplineX
#include "BSplineX/bspline/bspline_types.hpp"

namespace bsplinex::factory
{

template <typename T = double>
inline types::OpenUniform<T>
open_uniform(size_t degree, T begin, T end, size_t num_elems, std::vector<T> const &ctrl_points)
{
  return types::OpenUniform<T>{{begin, end, num_elems}, {ctrl_points}, degree};
}

template <typename T = double>
inline types::OpenUniform<T> open_uniform(size_t degree, T begin, T end, size_t num_elems)
{
  return open_uniform<T>(
      degree, begin, end, num_elems, std::vector<T>(num_elems - degree - 1, 0.0)
  );
}

template <typename T = double>
inline types::OpenUniform<T> open_uniform(size_t degree)
{
  T const begin{(T)0};
  T const end{(T)1};
  size_t const num_elems{1 + 2 * degree};
  return open_uniform<T>(
      degree, begin, end, num_elems, std::vector<T>(num_elems - degree - 1, 0.0)
  );
}

template <typename T = double>
inline types::OpenUniformConstant<T> open_uniform_constant(
    size_t degree, T begin, T end, size_t num_elems, std::vector<T> const &ctrl_points
)
{
  return types::OpenUniformConstant<T>{{begin, end, num_elems}, {ctrl_points}, degree};
}

template <typename T = double>
inline types::OpenUniformConstant<T>
open_uniform_constant(size_t degree, T begin, T end, size_t num_elems)
{
  return open_uniform_constant<T>(
      degree, begin, end, num_elems, std::vector<T>(num_elems - degree - 1, 0.0)
  );
}

template <typename T = double>
inline types::OpenUniformConstant<T> open_uniform_constant(size_t degree)
{
  T const begin{(T)0};
  T const end{(T)1};
  size_t const num_elems{1 + 2 * degree};
  return open_uniform_constant<T>(
      degree, begin, end, num_elems, std::vector<T>(num_elems - degree - 1, 0.0)
  );
}

template <typename T = double>
inline types::OpenNonUniform<T>
open_nonuniform(size_t degree, std::vector<T> const &knots, std::vector<T> const &ctrl_points)
{
  return types::OpenNonUniform<T>{{knots}, {ctrl_points}, degree};
}

template <typename T = double>
inline types::OpenNonUniform<T> open_nonuniform(size_t degree, std::vector<T> const &knots)
{
  return open_nonuniform<T>(degree, knots, std::vector<T>(knots.size() - degree - 1, 0.0));
}

template <typename T = double>
inline types::OpenNonUniform<T> open_nonuniform(size_t degree)
{
  std::vector<T> const knots(1 + 2 * degree);
  return open_nonuniform<T>(degree, knots, std::vector<T>(knots.size() - degree - 1, 0.0));
}

template <typename T = double>
inline types::OpenNonUniformConstant<T> open_nonuniform_constant(
    size_t degree, std::vector<T> const &knots, std::vector<T> const &ctrl_points
)
{
  return types::OpenNonUniformConstant<T>{{knots}, {ctrl_points}, degree};
}

template <typename T = double>
inline types::OpenNonUniformConstant<T>
open_nonuniform_constant(size_t degree, std::vector<T> const &knots)
{
  return open_nonuniform_constant<T>(degree, knots, std::vector<T>(knots.size() - degree - 1, 0.0));
}

template <typename T = double>
inline types::OpenNonUniformConstant<T> open_nonuniform_constant(size_t degree)
{
  std::vector<T> const knots(1 + 2 * degree);
  return open_nonuniform_constant<T>(degree, knots, std::vector<T>(knots.size() - degree - 1, 0.0));
}

} // namespace bsplinex::factory

#endif
