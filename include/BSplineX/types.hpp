#ifndef BSPLINEX_TYPES_HPP
#define BSPLINEX_TYPES_HPP

namespace bsplinex
{

/**
 * @brief Enum class representing different types of boundary conditions.
 */
enum class BoundaryCondition
{
  CLAMPED  = 0,
  OPEN     = 1,
  PERIODIC = 2
};

/**
 * @brief Enum class representing different types of curves.
 */
enum class Curve
{
  NON_UNIFORM = 0,
  UNIFORM     = 1
};

/**
 * @brief Enum class representing different types of extrapolation methods.
 */
enum class Extrapolation
{
  CONSTANT = 0,
  PERIODIC = 1,
  NONE     = 2
};

} // namespace bsplinex

#endif
