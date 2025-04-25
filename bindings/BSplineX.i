%module BSplineX

// Injected C++ code
%{
#include "BSplineX/bsplinex.hpp"
#include "BSplineX/defines.hpp"
#include "BSplineX/types.hpp"
#include "BSplineX/bspline/bspline.hpp"
#include "BSplineX/bspline/bspline_types.hpp"
#include "BSplineX/bspline/bspline_factory_open.hpp"
#include "BSplineX/bspline/bspline_factory_clamped.hpp"
#include "BSplineX/bspline/bspline_factory_periodic.hpp"
#include "BSplineX/bspline/bspline_lsq.hpp"
using namespace bsplinex;
%}

%include "std_pair.i"
%include "std_vector.i"
%template() std::pair<double, double>;
%template() std::vector<double>;

%include "BSplineX/bsplinex.hpp"
%include "BSplineX/defines.hpp"
%include "BSplineX/types.hpp"
%include "BSplineX/bspline/bspline.hpp"
%include "BSplineX/bspline/bspline_types.hpp"
%include "BSplineX/bspline/bspline_factory_open.hpp"
%include "BSplineX/bspline/bspline_factory_clamped.hpp"
%include "BSplineX/bspline/bspline_factory_periodic.hpp"
%include "BSplineX/bspline/bspline_lsq.hpp"
%include "exception.i"
using namespace bsplinex;

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%template(InterpolantCondition) InterpolantCondition<double>;
%template() std::vector<InterpolantCondition<double>>;

%template(OpenUniform) bspline::BSpline<double, Curve::UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE>;
%template() OpenUniform<double>;
%template(make_open_uniform) make_open_uniform<double>;

%template(OpenUniformConstant) bspline::BSpline<double, Curve::UNIFORM, BoundaryCondition::OPEN, Extrapolation::CONSTANT>;
%template() OpenUniformConstant<double>;
%template(make_open_uniform_constant) make_open_uniform_constant<double>;

%template(OpenNonUniform) bspline::BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE>;
%template() OpenNonUniform<double>;
%template(make_open_nonuniform) make_open_nonuniform<double>;

%template(OpenNonUniformConstant) bspline::BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::OPEN, Extrapolation::CONSTANT>;
%template() OpenNonUniformConstant<double>;
%template(make_open_nonuniform_constant) make_open_nonuniform_constant<double>;

%template(ClampedUniform) bspline::BSpline<double, Curve::UNIFORM, BoundaryCondition::CLAMPED, Extrapolation::NONE>;
%template() ClampedUniform<double>;
%template(make_clamped_uniform) make_clamped_uniform<double>;

%template(ClampedUniformConstant) bspline::BSpline<double, Curve::UNIFORM, BoundaryCondition::CLAMPED, Extrapolation::CONSTANT>;
%template() ClampedUniformConstant<double>;
%template(make_clamped_uniform_constant) make_clamped_uniform_constant<double>;

%template(ClampedNonUniform) bspline::BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::CLAMPED, Extrapolation::NONE>;
%template() ClampedNonUniform<double>;
%template(make_clamped_nonuniform) make_clamped_nonuniform<double>;

%template(ClampedNonUniformConstant) bspline::BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::CLAMPED, Extrapolation::CONSTANT>;
%template() ClampedNonUniformConstant<double>;
%template(make_clamped_nonuniform_constant) make_clamped_nonuniform_constant<double>;

%template(PeriodicUniform) bspline::BSpline<double, Curve::UNIFORM, BoundaryCondition::PERIODIC, Extrapolation::PERIODIC>;
%template() PeriodicUniform<double>;
%template(make_periodic_uniform) make_periodic_uniform<double>;

%template(PeriodicNonUniform) bspline::BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::PERIODIC, Extrapolation::PERIODIC>;
%template() PeriodicNonUniform<double>;
%template(make_periodic_nonuniform) make_periodic_nonuniform<double>;
