#ifndef BSPLINEX_TYPES_HPP
#define BSPLINEX_TYPES_HPP

namespace bsplinex
{

/**
 * @brief Enum class representing different types of boundary conditions.
 */
enum class BoundaryCondition : uint8_t
{
  CLAMPED  = 0,
  OPEN     = 1,
  PERIODIC = 2
};

/**
 * @brief Enum class representing different types of curves.
 */
enum class Curve : uint8_t
{
  NON_UNIFORM = 0,
  UNIFORM     = 1
};

/**
 * @brief Enum class representing different types of extrapolation methods.
 */
enum class Extrapolation : uint8_t
{
  CONSTANT = 0,
  PERIODIC = 1,
  NONE     = 2
};

} // namespace bsplinex

#endif
