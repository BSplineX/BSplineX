#ifndef BSPLINEX_BSPLINE_BSPLINE_FACTORY_PERIODIC_HPP
#define BSPLINEX_BSPLINE_BSPLINE_FACTORY_PERIODIC_HPP

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
 * @brief Creates a periodic uniform BSpline.
 *
 * This is a generic constructor in case you know exactly the parameters of the
 * BSpline.
 * Note that this being a periodic BSpline, the knot and control points vectors
 * you specify will be respected, but some padding will be applied to ensure
 * periodicity.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @param begin The begin of the uniform knot vector.
 * @param end The end of the uniform knot vector.
 * @param num_elems The number of elements in the knot vector.
 * @param ctrl_points The control points of the BSpline.
 * @return A periodic uniform BSpline.
 */
template <typename T = double>
types::PeriodicUniform<T> make_periodic_uniform(
    size_t degree, T begin, T end, size_t num_elems, std::vector<T> const &ctrl_points
)
{
  return types::PeriodicUniform<T>{
      knots::Data<T, types::PeriodicUniform<T>::curve_type>{begin, end, num_elems},
      control_points::Data<T>{ctrl_points},
      degree
  };
}

/**
 * @brief Creates a periodic uniform BSpline.
 *
 * This is a useful constructor in case you want to fit a BSpline to some data.
 * Note that this being a periodic BSpline, the knot vector you specify will be
 * respected, but some padding will be applied to ensure periodicity.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @param begin The begin of the uniform knot vector.
 * @param end The end of the uniform knot vector.
 * @param num_elems The number of elements in the knot vector.
 * @return A periodic uniform BSpline.
 */
template <typename T = double>
types::PeriodicUniform<T> make_periodic_uniform(size_t degree, T begin, T end, size_t num_elems)
{
  return make_periodic_uniform<T>(degree, begin, end, num_elems, std::vector<T>(num_elems - 1));
}

/**
 * @brief Creates a periodic uniform BSpline.
 *
 * This is a useful constructor in case you want to interpolate a BSpline to
 * some data.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @return A periodic uniform BSpline.
 */
template <typename T = double>
types::PeriodicUniform<T> make_periodic_uniform(size_t degree)
{
  constexpr T begin{ZERO<T>};
  constexpr T end{ONE<T>};
  size_t const num_elems{2 * degree};
  return make_periodic_uniform<T>(degree, begin, end, num_elems, std::vector<T>(num_elems - 1));
}

/**
 * @brief Creates a periodic non-uniform BSpline.
 *
 * This is a generic constructor in case you know exactly the parameters of the
 * BSpline.
 * Note that this being a periodic BSpline, the knot and control points vectors
 * you specify will be respected, but some padding will be applied to ensure
 * periodicity.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @param knots The knot vector of the BSpline.
 * @param ctrl_points The control points of the BSpline.
 * @return A periodic non-uniform BSpline.
 */
template <typename T = double>
types::PeriodicNonUniform<T> make_periodic_nonuniform(
    size_t degree, std::vector<T> const &knots, std::vector<T> const &ctrl_points
)
{
  return types::PeriodicNonUniform<T>{
      knots::Data<T, types::PeriodicNonUniform<T>::curve_type>{knots},
      control_points::Data<T>{ctrl_points},
      degree
  };
}

/**
 * @brief Creates a periodic non-uniform BSpline.
 *
 * This is a useful constructor in case you want to fit a BSpline to some data.
 * Note that this being a periodic BSpline, the knot vector you specify will be
 * respected, but some padding will be applied to ensure periodicity.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @param knots The knot vector of the BSpline.
 * @return A periodic non-uniform BSpline.
 */
template <typename T = double>
types::PeriodicNonUniform<T> make_periodic_nonuniform(size_t degree, std::vector<T> const &knots)
{
  return make_periodic_nonuniform<T>(degree, knots, std::vector<T>(knots.size() - 1));
}

/**
 * @brief Creates a periodic non-uniform BSpline.
 *
 * This is a useful constructor in case you want to interpolate a BSpline to
 * some data.
 *
 * @tparam T The numeric type.
 * @param degree The degree of the BSpline.
 * @return A periodic non-uniform BSpline.
 */
template <typename T = double>
types::PeriodicNonUniform<T> make_periodic_nonuniform(size_t degree)
{
  std::vector<T> const knots(2 * degree);
  return make_periodic_nonuniform<T>(degree, knots, std::vector<T>(knots.size() - 1));
}

} // namespace bsplinex::factory

#endif
