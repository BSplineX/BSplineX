#ifndef BSPLINEX_KNOTS_KNOTS_HPP
#define BSPLINEX_KNOTS_KNOTS_HPP

// Standard includes
#include <utility>

// BSplineX includes
#include "BSplineX/knots/t_atter.hpp"
#include "BSplineX/knots/t_extrapolator.hpp"
#include "BSplineX/knots/t_finder.hpp"
#include "BSplineX/types.hpp"
#include "BSplineX/windows.hpp"

/**
 * Naming convention:
 * - `m` -> number of knots
 * - `p` -> degree of the curve
 * - `t` -> knots vector
 *
 * Curve domain:
 * - If the curve is open, the domain is [t_p, t_{end - p}]
 * - If the curve is periodic, the domain is [t_0, t_{end}] but appropriate
 *   padding is needed
 * - If the curve is clamped, the domain is [t_0, t_{end}] but the start and end
 *   knots must have multiplicity `p+1`
 *
 * Knots padding:
 * - If the curve is open, no padding is needed, the full `n + p + 1` knots have
 *   to be provided
 * - If the curve is periodic, we need to add `p` knots at the left and right
 *   following periodicity: [0, 1, 2, 2.5, 3] with p = 3 ->
 *   [-2.0, -1.0, -0.5, 0, 1, 2, 2.5, 3, 4, 5, 5.5]
 * - If the curve is clamped, we must repeat the first and last knots `p` times:
 *   [0, 1, 2, 2.5, 3] with p = 3 -> [0, 0, 0, 0, 1, 2, 2.5, 3, 3, 3, 3]
 *
 */

namespace bsplinex::knots
{

template <typename T, Curve C, BoundaryCondition BC, Extrapolation EXT>
class Knots
{
private:
  Atter<T, C, BC> m_atter{};
  Extrapolator<T, C, BC, EXT> m_extrapolator{};
  T m_value_left{};
  T m_value_right{};
  size_t m_degree{};

public:
  Knots() = default;

  Knots(Data<T, C> const &data, size_t degree)
      : m_atter{data, degree}, m_extrapolator{this->m_atter, degree},
        m_value_left{this->m_atter.at(degree)},
        m_value_right{this->m_atter.at(this->m_atter.size() - degree - 1)}, m_degree{degree}
  {
  }

  Knots(Knots const &other) = default;

  Knots(Knots &&other) = default;

  ~Knots() = default;

  Knots &operator=(Knots const &other) = default;

  Knots &operator=(Knots &&other) = default;

  [[nodiscard]] std::pair<size_t, T> find(T value) const
  {
    if (value < this->m_value_left or value > this->m_value_right)
    {
      value = this->m_extrapolator.extrapolate(value);
    }

    return std::pair<size_t, T>{
        knots::find<T, C, BC, EXT>(this->m_atter, this->m_degree, value), value
    };
  }

  [[nodiscard]] std::pair<T, T> domain() const
  {
    return std::make_pair(m_value_left, m_value_right);
  }

  [[nodiscard]] T at(size_t index) const { return this->m_atter.at(index); }

  [[nodiscard]] size_t size() const { return this->m_atter.size(); }

  [[nodiscard]] Knots get_derivative_knots() const
  {
    Atter<T, C, BC> d_atter = this->m_atter;
    return Knots(d_atter.pop_tails(), this->m_degree - 1);
  }

private:
  Knots(Atter<T, C, BC> const &atter, size_t degree)
      : m_atter{atter}, m_extrapolator{this->m_atter, degree},
        m_value_left{this->m_atter.at(degree)},
        m_value_right{this->m_atter.at(this->m_atter.size() - degree - 1)}, m_degree{degree}
  {
  }
};

} // namespace bsplinex::knots

#endif
