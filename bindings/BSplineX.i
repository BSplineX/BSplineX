%module BSplineX

// Injected C++ code
%{
#include "BSplineX/bspline/bspline.hpp"
#include "BSplineX/bspline/bspline_types.hpp"
#include "BSplineX/bspline/bspline_factory.hpp"
using namespace bsplinex;
%}

%include "std_pair.i"
%include "std_vector.i"
%template() std::pair<double, double>;
%template() std::vector<double>;


%include "BSplineX/bspline/bspline.hpp"
%include "BSplineX/bspline/bspline_factory.hpp"
%include "BSplineX/bspline/bspline_types.hpp"
%include "exception.i"
using namespace bsplinex;

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

//  ██████╗ ██████╗ ███████╗███╗   ██╗
// ██╔═══██╗██╔══██╗██╔════╝████╗  ██║
// ██║   ██║██████╔╝█████╗  ██╔██╗ ██║
// ██║   ██║██╔═══╝ ██╔══╝  ██║╚██╗██║
// ╚██████╔╝██║     ███████╗██║ ╚████║
//  ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═══╝

%template(OpenUniform) bspline::BSpline<double, Curve::UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE>;
%template() types::OpenUniform<double>;
%template(open_uniform) factory::open_uniform<double>;

%template(OpenUniformConstant) bspline::BSpline<double, Curve::UNIFORM, BoundaryCondition::OPEN, Extrapolation::CONSTANT>;
%template() types::OpenUniformConstant<double>;
%template(open_uniform_constant) factory::open_uniform_constant<double>;

%template(OpenNonUniform) bspline::BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE>;
%template() types::OpenNonUniform<double>;
%template(open_nonuniform) factory::open_nonuniform<double>;

%template(OpenNonUniformConstant) bspline::BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::OPEN, Extrapolation::CONSTANT>;
%template() types::OpenNonUniformConstant<double>;
%template(open_nonuniform_constant) factory::open_nonuniform_constant<double>;

//  ██████╗██╗      █████╗ ███╗   ███╗██████╗ ███████╗██████╗
// ██╔════╝██║     ██╔══██╗████╗ ████║██╔══██╗██╔════╝██╔══██╗
// ██║     ██║     ███████║██╔████╔██║██████╔╝█████╗  ██║  ██║
// ██║     ██║     ██╔══██║██║╚██╔╝██║██╔═══╝ ██╔══╝  ██║  ██║
// ╚██████╗███████╗██║  ██║██║ ╚═╝ ██║██║     ███████╗██████╔╝
// ╚═════╝╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝     ╚══════╝╚═════╝

%template(ClampedUniform) bspline::BSpline<double, Curve::UNIFORM, BoundaryCondition::CLAMPED, Extrapolation::NONE>;
%template() types::ClampedUniform<double>;
%template(clamped_uniform) factory::clamped_uniform<double>;

%template(ClampedUniformConstant) bspline::BSpline<double, Curve::UNIFORM, BoundaryCondition::CLAMPED, Extrapolation::CONSTANT>;
%template() types::ClampedUniformConstant<double>;
%template(clamped_uniform_constant) factory::clamped_uniform_constant<double>;

%template(ClampedNonUniform) bspline::BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::CLAMPED, Extrapolation::NONE>;
%template() types::ClampedNonUniform<double>;
%template(clamped_nonuniform) factory::clamped_nonuniform<double>;

%template(ClampedNonUniformConstant) bspline::BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::CLAMPED, Extrapolation::CONSTANT>;
%template() types::ClampedNonUniformConstant<double>;
%template(clamped_nonuniform_constant) factory::clamped_nonuniform_constant<double>;

// ██████╗ ███████╗██████╗ ██╗ ██████╗ ██████╗ ██╗ ██████╗
// ██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██╔══██╗██║██╔════╝
// ██████╔╝█████╗  ██████╔╝██║██║   ██║██║  ██║██║██║
// ██╔═══╝ ██╔══╝  ██╔══██╗██║██║   ██║██║  ██║██║██║
// ██║     ███████╗██║  ██║██║╚██████╔╝██████╔╝██║╚██████╗
// ╚═╝     ╚══════╝╚═╝  ╚═╝╚═╝ ╚═════╝ ╚═════╝ ╚═╝ ╚═════╝

%template(PeriodicUniform) bspline::BSpline<double, Curve::UNIFORM, BoundaryCondition::PERIODIC, Extrapolation::PERIODIC>;
%template() types::PeriodicUniform<double>;
%template(periodic_uniform) factory::periodic_uniform<double>;

%template(PeriodicNonUniform) bspline::BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::PERIODIC, Extrapolation::PERIODIC>;
%template() types::PeriodicNonUniform<double>;
%template(periodic_nonuniform) factory::periodic_nonuniform<double>;

