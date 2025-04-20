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

using factory::make_open_nonuniform;
using factory::make_open_nonuniform_constant;
using factory::make_open_uniform;
using factory::make_open_uniform_constant;

using types::ClampedNonUniform;
using types::ClampedNonUniformConstant;
using types::ClampedUniform;
using types::ClampedUniformConstant;

using factory::make_clamped_nonuniform;
using factory::make_clamped_nonuniform_constant;
using factory::make_clamped_uniform;
using factory::make_clamped_uniform_constant;

using types::PeriodicNonUniform;
using types::PeriodicUniform;

using factory::make_periodic_nonuniform;
using factory::make_periodic_uniform;

} // namespace bsplinex

#endif
