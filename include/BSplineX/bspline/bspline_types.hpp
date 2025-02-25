#ifndef BSPLINEX_BSPLINE_BSPLINE_TYPES_HPP
#define BSPLINEX_BSPLINE_BSPLINE_TYPES_HPP

#include "BSplineX/bspline/bspline.hpp"

namespace bsplinex::types
{

template <typename T = double>
using OpenUniform =
    bspline::BSpline<T, Curve::UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE>;

template <typename T = double>
using OpenUniformConstant =
    bspline::BSpline<T, Curve::UNIFORM, BoundaryCondition::OPEN, Extrapolation::CONSTANT>;

template <typename T = double>
using OpenNonUniform =
    bspline::BSpline<T, Curve::NON_UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE>;

template <typename T = double>
using OpenNonUniformConstant =
    bspline::BSpline<T, Curve::NON_UNIFORM, BoundaryCondition::OPEN, Extrapolation::CONSTANT>;

template <typename T = double>
using ClampedUniform =
    bspline::BSpline<T, Curve::UNIFORM, BoundaryCondition::CLAMPED, Extrapolation::NONE>;

template <typename T = double>
using ClampedUniformConstant =
    bspline::BSpline<T, Curve::UNIFORM, BoundaryCondition::CLAMPED, Extrapolation::CONSTANT>;

template <typename T = double>
using ClampedNonUniform =
    bspline::BSpline<T, Curve::NON_UNIFORM, BoundaryCondition::CLAMPED, Extrapolation::NONE>;

template <typename T = double>
using ClampedNonUniformConstant =
    bspline::BSpline<T, Curve::NON_UNIFORM, BoundaryCondition::CLAMPED, Extrapolation::CONSTANT>;

template <typename T = double>
using PeriodicUniform =
    bspline::BSpline<T, Curve::UNIFORM, BoundaryCondition::PERIODIC, Extrapolation::PERIODIC>;

template <typename T = double>
using PeriodicNonUniform =
    bspline::BSpline<T, Curve::NON_UNIFORM, BoundaryCondition::PERIODIC, Extrapolation::PERIODIC>;

} // namespace bsplinex::types

#endif
