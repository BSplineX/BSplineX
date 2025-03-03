#ifndef BSPLINEX_CONTROL_POINTS_CONTROL_POINTS_HPP
#define BSPLINEX_CONTROL_POINTS_CONTROL_POINTS_HPP

// Standard includes
#include <cstddef>

// BSplineX includes
#include "BSplineX/control_points/c_atter.hpp"
#include "BSplineX/defines.hpp"
#include "BSplineX/types.hpp"

/**
 * Naming convention:
 * - `n` -> number on control points
 * - `m` -> number of knots
 * - `p` -> degree of the curve
 * - `c` -> control points vector
 *
 * Control points:
 * - If the curve is open, the number of control points must be `n = m - p - 1`
 * - If the curve is periodic, the number of control points must be `n = m - 1`
 *   since they will need to be padded for periodicity by repeating the first
 *   `p` control points at the end
 * - If the curve is clamped, the number of control points must be
 *   `n = m + p - 1`
 *
 */

namespace bsplinex::control_points
{

template <typename T, BoundaryCondition BC>
class ControlPoints
{
private:
  Atter<T, BC> atter{};
  size_t degree{};

public:
  ControlPoints() { DEBUG_LOG_CALL(); }

  ControlPoints(Data<T> data, size_t degree) : atter{data, degree}, degree{degree}
  {
    DEBUG_LOG_CALL();
  }

  ControlPoints(ControlPoints const &other) : atter(other.atter), degree(other.degree)
  {
    DEBUG_LOG_CALL();
  }

  ControlPoints(ControlPoints &&other) noexcept
      : atter(std::move(other.atter)), degree(other.degree)
  {
    DEBUG_LOG_CALL();
  }

  ~ControlPoints() noexcept { DEBUG_LOG_CALL(); }

  ControlPoints &operator=(ControlPoints const &other)
  {
    DEBUG_LOG_CALL();
    if (this == &other)
      return *this;
    atter  = other.atter;
    degree = other.degree;
    return *this;
  }

  ControlPoints &operator=(ControlPoints &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
      return *this;
    atter  = std::move(other.atter);
    degree = other.degree;
    return *this;
  }

  T at(size_t index) const { return this->atter.at(index); }

  [[nodiscard]] size_t size() const { return this->atter.size(); }
  
};

} // namespace bsplinex::control_points

#endif