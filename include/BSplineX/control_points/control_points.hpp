#ifndef BSPLINEX_CONTROL_POINTS_CONTROL_POINTS_HPP
#define BSPLINEX_CONTROL_POINTS_CONTROL_POINTS_HPP

// Standard includes
#include <cstddef>

// BSplineX includes
#include "BSplineX/control_points/c_atter.hpp"
#include "BSplineX/knots/knots.hpp"
#include "BSplineX/types.hpp"
#include "BSplineX/windows.hpp"

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
  ControlPoints() = default;

  ControlPoints(Data<T> data, size_t degree) : atter{data, degree}, degree{degree} {}

  ControlPoints(ControlPoints const &other) = default;

  ControlPoints(ControlPoints &&other) noexcept = default;

  ~ControlPoints() noexcept = default;

  ControlPoints &operator=(ControlPoints const &other) = default;

  ControlPoints &operator=(ControlPoints &&other) noexcept = default;

  [[nodiscard]] T at(size_t index) const { return this->atter.at(index); }

  [[nodiscard]] size_t size() const { return this->atter.size(); }

  template <Curve C, Extrapolation EXT>
  [[nodiscard]] ControlPoints
  get_derivative_control_points(knots::Knots<T, C, BC, EXT> const &knots) const
  {
    size_t const d_num_ctrl_pts = this->atter.get_derivative_data_size();
    std::vector<T> d_ctrl_points;
    d_ctrl_points.reserve(d_num_ctrl_pts);
    for (size_t i = 0; i < d_num_ctrl_pts; i++)
    {
      d_ctrl_points.push_back(
          static_cast<T>(this->degree) / (knots.at(i + this->degree + 1) - knots.at(i + 1)) *
          (this->at(i + 1) - this->at(i))
      );
    }

    return ControlPoints({d_ctrl_points}, this->degree - 1);
  }
};

} // namespace bsplinex::control_points

#endif
