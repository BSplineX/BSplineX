#ifndef BSPLINE_FACTORY_CLAMPED_HPP
#define BSPLINE_FACTORY_CLAMPED_HPP

// Standard
#include <vector>

// BSplineX
#include "BSplineX/bspline/bspline_types.hpp"

namespace bsplinex::factory
{

template <typename T = double>
inline types::ClampedUniform<T>
clamped_uniform(size_t degree, T begin, T end, size_t num_elems, std::vector<T> const &ctrl_points)
{
  return types::ClampedUniform<T>{{begin, end, num_elems}, {ctrl_points}, degree};
}

template <typename T = double>
inline types::ClampedUniform<T> clamped_uniform(size_t degree, T begin, T end, size_t num_elems)
{
  return clamped_uniform<T>(
      degree, begin, end, num_elems, std::vector<T>(num_elems + degree - 1, 0.0)
  );
}

template <typename T = double>
inline types::ClampedUniform<T> clamped_uniform(size_t degree)
{
  T const begin{(T)0};
  T const end{(T)1};
  size_t const num_elems{2};
  return clamped_uniform<T>(
      degree, begin, end, num_elems, std::vector<T>(num_elems + degree - 1, 0.0)
  );
}

template <typename T = double>
inline types::ClampedUniformConstant<T> clamped_uniform_constant(
    size_t degree, T begin, T end, size_t num_elems, std::vector<T> const &ctrl_points
)
{
  return types::ClampedUniformConstant<T>{{begin, end, num_elems}, {ctrl_points}, degree};
}

template <typename T = double>
inline types::ClampedUniformConstant<T>
clamped_uniform_constant(size_t degree, T begin, T end, size_t num_elems)
{
  return clamped_uniform_constant<T>(
      degree, begin, end, num_elems, std::vector<T>(num_elems + degree - 1, 0.0)
  );
}

template <typename T = double>
inline types::ClampedUniformConstant<T> clamped_uniform_constant(size_t degree)
{
  T const begin{(T)0};
  T const end{(T)1};
  size_t const num_elems{2};
  return clamped_uniform_constant<T>(
      degree, begin, end, num_elems, std::vector<T>(num_elems + degree - 1, 0.0)
  );
}

template <typename T = double>
inline types::ClampedNonUniform<T>
clamped_nonuniform(size_t degree, std::vector<T> const &knots, std::vector<T> const &ctrl_points)
{
  return types::ClampedNonUniform<T>{{knots}, {ctrl_points}, degree};
}

template <typename T = double>
inline types::ClampedNonUniform<T> clamped_nonuniform(size_t degree, std::vector<T> const &knots)
{
  return clamped_nonuniform<T>(degree, knots, std::vector<T>(knots.size() + degree - 1, 0.0));
}

template <typename T = double>
inline types::ClampedNonUniform<T> clamped_nonuniform(size_t degree)
{
  std::vector<T> const knots(2);
  return clamped_nonuniform<T>(degree, knots, std::vector<T>(knots.size() + degree - 1, 0.0));
}

template <typename T = double>
inline types::ClampedNonUniformConstant<T> clamped_nonuniform_constant(
    size_t degree, std::vector<T> const &knots, std::vector<T> const &ctrl_points
)
{
  return types::ClampedNonUniformConstant<T>{{knots}, {ctrl_points}, degree};
}

template <typename T = double>
inline types::ClampedNonUniformConstant<T>
clamped_nonuniform_constant(size_t degree, std::vector<T> const &knots)
{
  return clamped_nonuniform_constant<T>(
      degree, knots, std::vector<T>(knots.size() + degree - 1, 0.0)
  );
}

template <typename T = double>
inline types::ClampedNonUniformConstant<T> clamped_nonuniform_constant(size_t degree)
{
  std::vector<T> const knots(2);
  return clamped_nonuniform_constant<T>(
      degree, knots, std::vector<T>(knots.size() + degree - 1, 0.0)
  );
}

} // namespace bsplinex::factory

#endif
