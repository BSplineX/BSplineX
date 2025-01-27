#ifndef BSPLINEX_TYPES_HPP
#define BSPLINEX_TYPES_HPP

namespace bsplinex
{

enum class BoundaryCondition
{
  CLAMPED  = 0,
  OPEN     = 1,
  PERIODIC = 2
};

enum class Curve
{
  NON_UNIFORM = 0,
  UNIFORM     = 1
};

enum class Extrapolation
{
  CONSTANT = 0,
  PERIODIC = 1,
  NONE     = 2
};

} // namespace bsplinex

#endif
