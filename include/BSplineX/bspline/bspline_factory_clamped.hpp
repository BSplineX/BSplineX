#ifndef BSPLINEX_BSPLINE_BSPLINE_FACTORY_CLAMPED_HPP
#define BSPLINEX_BSPLINE_BSPLINE_FACTORY_CLAMPED_HPP

// Standard
#include <vector>

// BSplineX
#include "BSplineX/bspline/bspline_types.hpp"
#include "BSplineX/control_points/c_data.hpp"
#include "BSplineX/defines.hpp"
#include "BSplineX/knots/t_data.hpp"

namespace bsplinex::factory
{

using namespace constants;

/**
 * @brief Creates a clamped uniform BSpline.
 *
 * This is a generic constructor in case you know exactly the parameters of the
 * BSpline.
 * Note that this being a clamped BSpline, the knot vector you specify will be
 * respected and the first and last knots will have multiplicity `p+1`.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @param begin The begin of the uniform knot vector.
 * @param end The end of the uniform knot vector.
 * @param num_elems The number of elements in the knot vector.
 * @param ctrl_points The control points of the BSpline.
 * @return A clamped uniform BSpline.
 */
template <typename T = double>
inline types::ClampedUniform<T>
clamped_uniform(size_t degree, T begin, T end, size_t num_elems, std::vector<T> const &ctrl_points)
{
  return types::ClampedUniform<T>{
      knots::Data<T, types::ClampedUniform<T>::curve_type>{begin, end, num_elems},
      control_points::Data<T>{ctrl_points},
      degree
  };
}

/**
 * @brief Creates a clamped uniform BSpline.
 *
 * This is a useful constructor in case you want to fit a BSpline to some data.
 * Note that this being a clamped BSpline, the knot vector you specify will be
 * respected and the first and last knots will have multiplicity `p+1`.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @param begin The begin of the uniform knot vector.
 * @param end The end of the uniform knot vector.
 * @param num_elems The number of elements in the knot vector.
 * @return A clamped uniform BSpline.
 */
template <typename T = double>
inline types::ClampedUniform<T> clamped_uniform(size_t degree, T begin, T end, size_t num_elems)
{
  return clamped_uniform<T>(degree, begin, end, num_elems, std::vector<T>(num_elems + degree - 1));
}

/**
 * @brief Creates a clamped uniform BSpline.
 *
 * This is a useful constructor in case you want to interpolate a BSpline to
 * some data.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @return A clamped uniform BSpline.
 */
template <typename T = double>
inline types::ClampedUniform<T> clamped_uniform(size_t degree)
{
  constexpr T begin{ZERO<T>};
  constexpr T end{ONE<T>};
  constexpr size_t num_elems{2};
  return clamped_uniform<T>(degree, begin, end, num_elems, std::vector<T>(num_elems + degree - 1));
}

/**
 * @brief Creates a clamped uniform BSpline with constant extrapolation.
 *
 * This is a generic constructor in case you know exactly the parameters of the
 * BSpline.
 * Note that this being a clamped BSpline, the knot vector you specify will be
 * respected and the first and last knots will have multiplicity `p+1`.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @param begin The begin of the uniform knot vector.
 * @param end The end of the uniform knot vector.
 * @param num_elems The number of elements in the knot vector.
 * @param ctrl_points The control points of the BSpline.
 * @return A clamped uniform constant BSpline.
 */
template <typename T = double>
inline types::ClampedUniformConstant<T> clamped_uniform_constant(
    size_t degree, T begin, T end, size_t num_elems, std::vector<T> const &ctrl_points
)
{
  return types::ClampedUniformConstant<T>{
      knots::Data<T, types::ClampedUniform<T>::curve_type>{begin, end, num_elems},
      control_points::Data<T>{ctrl_points},
      degree
  };
}

/**
 * @brief Creates a clamped uniform BSpline with constant extrapolation.
 *
 * This is a useful constructor in case you want to fit a BSpline to some data.
 * Note that this being a clamped BSpline, the knot vector you specify will be
 * respected and the first and last knots will have multiplicity `p+1`.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @param begin The begin of the uniform knot vector.
 * @param end The end of the uniform knot vector.
 * @param num_elems The number of elements in the knot vector.
 * @return A clamped uniform constant BSpline.
 */
template <typename T = double>
inline types::ClampedUniformConstant<T>
clamped_uniform_constant(size_t degree, T begin, T end, size_t num_elems)
{
  return clamped_uniform_constant<T>(
      degree, begin, end, num_elems, std::vector<T>(num_elems + degree - 1)
  );
}

/**
 * @brief Creates a clamped uniform BSpline with constant extrapolation.
 *
 * This is a useful constructor in case you want to interpolate a BSpline to
 * some data.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @return A clamped uniform constant BSpline.
 */
template <typename T = double>
inline types::ClampedUniformConstant<T> clamped_uniform_constant(size_t degree)
{
  constexpr T begin{ZERO<T>};
  constexpr T end{ONE<T>};
  constexpr size_t num_elems{2};
  return clamped_uniform_constant<T>(
      degree, begin, end, num_elems, std::vector<T>(num_elems + degree - 1)
  );
}

/**
 * @brief Creates a clamped non-uniform BSpline.
 *
 * This is a generic constructor in case you know exactly the parameters of the
 * BSpline.
 * Note that this being a clamped BSpline, the knot vector you specify will be
 * respected and the first and last knots will have multiplicity `p+1`.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @param knots The knot vector of the BSpline.
 * @param ctrl_points The control points of the BSpline.
 * @return A clamped non-uniform BSpline.
 */
template <typename T = double>
inline types::ClampedNonUniform<T>
clamped_nonuniform(size_t degree, std::vector<T> const &knots, std::vector<T> const &ctrl_points)
{
  return types::ClampedNonUniform<T>{
      knots::Data<T, types::ClampedNonUniform<T>::curve_type>{knots},
      control_points::Data<T>{ctrl_points},
      degree
  };
}

/**
 * @brief Creates a clamped non-uniform BSpline.
 *
 * This is a useful constructor in case you want to fit a BSpline to some data.
 * Note that this being a clamped BSpline, the knot vector you specify will be
 * respected and the first and last knots will have multiplicity `p+1`.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @param knots The knot vector of the BSpline.
 * @return A clamped non-uniform BSpline.
 */
template <typename T = double>
inline types::ClampedNonUniform<T> clamped_nonuniform(size_t degree, std::vector<T> const &knots)
{
  return clamped_nonuniform<T>(degree, knots, std::vector<T>(knots.size() + degree - 1));
}

/**
 * @brief Creates a clamped non-uniform BSpline.
 *
 * This is a useful constructor in case you want to interpolate a BSpline to
 * some data.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @return A clamped non-uniform BSpline.
 */
template <typename T = double>
inline types::ClampedNonUniform<T> clamped_nonuniform(size_t degree)
{
  std::vector<T> const knots(2);
  return clamped_nonuniform<T>(degree, knots, std::vector<T>(knots.size() + degree - 1));
}

/**
 * @brief Creates a clamped non-uniform BSpline with constant extrapolation.
 *
 * This is a generic constructor in case you know exactly the parameters of the
 * BSpline.
 * Note that this being a clamped BSpline, the knot vector you specify will be
 * respected and the first and last knots will have multiplicity `p+1`.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @param knots The knot vector of the BSpline.
 * @param ctrl_points The control points of the BSpline.
 * @return A clamped non-uniform constant BSpline.
 */
template <typename T = double>
inline types::ClampedNonUniformConstant<T> clamped_nonuniform_constant(
    size_t degree, std::vector<T> const &knots, std::vector<T> const &ctrl_points
)
{
  return types::ClampedNonUniformConstant<T>{
      knots::Data<T, types::ClampedNonUniform<T>::curve_type>{knots},
      control_points::Data<T>{ctrl_points},
      degree
  };
}

/**
 * @brief Creates a clamped non-uniform BSpline with constant extrapolation.
 *
 * This is a useful constructor in case you want to fit a BSpline to some data.
 * Note that this being a clamped BSpline, the knot vector you specify will be
 * respected and the first and last knots will have multiplicity `p+1`.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @param knots The knot vector of the BSpline.
 * @return A clamped non-uniform constant BSpline.
 */
template <typename T = double>
inline types::ClampedNonUniformConstant<T>
clamped_nonuniform_constant(size_t degree, std::vector<T> const &knots)
{
  return clamped_nonuniform_constant<T>(degree, knots, std::vector<T>(knots.size() + degree - 1));
}

/**
 * @brief Creates a clamped non-uniform BSpline with constant extrapolation.
 *
 * This is a useful constructor in case you want to interpolate a BSpline to
 * some data.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @return A clamped non-uniform constant BSpline.
 */
template <typename T = double>
inline types::ClampedNonUniformConstant<T> clamped_nonuniform_constant(size_t degree)
{
  std::vector<T> const knots(2);
  return clamped_nonuniform_constant<T>(degree, knots, std::vector<T>(knots.size() + degree - 1));
}

} // namespace bsplinex::factory

#endif
