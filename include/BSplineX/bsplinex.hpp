#ifndef BSPLINEX_BSPLINEX_HPP
#define BSPLINEX_BSPLINEX_HPP

#include "BSplineX/bspline/bspline_factory_clamped.hpp"
#include "BSplineX/bspline/bspline_factory_open.hpp"
#include "BSplineX/bspline/bspline_factory_periodic.hpp"
#include "BSplineX/bspline/bspline_types.hpp"

namespace bsplinex
{

using types::OpenNonUniform;
using types::OpenNonUniformConstant;
using types::OpenUniform;
using types::OpenUniformConstant;

using factory::open_nonuniform;
using factory::open_nonuniform_constant;
using factory::open_uniform;
using factory::open_uniform_constant;

using types::ClampedNonUniform;
using types::ClampedNonUniformConstant;
using types::ClampedUniform;
using types::ClampedUniformConstant;

using factory::clamped_nonuniform;
using factory::clamped_nonuniform_constant;
using factory::clamped_uniform;
using factory::clamped_uniform_constant;

using types::PeriodicNonUniform;
using types::PeriodicUniform;

using factory::periodic_nonuniform;
using factory::periodic_uniform;

} // namespace bsplinex

#endif
