#ifndef BSPLINEX_KNOTS_T_EXTRAPOLATOR_HPP
#define BSPLINEX_KNOTS_T_EXTRAPOLATOR_HPP

// Standard includes
#include <cmath>
#include <cstddef>

// BSplineX includes
#include "BSplineX/defines.hpp"
#include "BSplineX/knots/t_atter.hpp"
#include "BSplineX/types.hpp"
#include "BSplineX/windows.hpp"

namespace bsplinex::knots
{

template <typename T, Curve C, BoundaryCondition BC, Extrapolation EXT>
class Extrapolator
{
public:
  [[nodiscard]] virtual size_t extrapolate(T value) const = 0;
};

template <typename T, Curve C, BoundaryCondition BC>
class Extrapolator<T, C, BC, Extrapolation::NONE>
{
public:
  Extrapolator() = default;

  Extrapolator(Atter<T, C, BC> const & /*atter*/, size_t /*degree*/) {}

  [[nodiscard]] T extrapolate(T /*value*/) const
  {
    releaseassert(false, "Extrapolation explicitly set to NONE");

    return constants::ZERO<T>;
  }
};

template <typename T, Curve C, BoundaryCondition BC>
class Extrapolator<T, C, BC, Extrapolation::CONSTANT>
{
private:
  T m_value_left{};
  T m_value_right{};

public:
  Extrapolator() = default;

  Extrapolator(Atter<T, C, BC> const &atter, size_t degree)
      : m_value_left{atter.at(degree)}, m_value_right{atter.at(atter.size() - degree - 1)}
  {
  }

  Extrapolator(Extrapolator const &other) = default;

  Extrapolator(Extrapolator &&other) = default;

  ~Extrapolator() = default;

  Extrapolator &operator=(Extrapolator const &other) = default;

  Extrapolator &operator=(Extrapolator &&other) = default;

  [[nodiscard]] T extrapolate(T value) const
  {
    debugassert(
        value < this->m_value_left or value > this->m_value_right, "Value not outside of the domain"
    );
    return value < this->m_value_left ? this->m_value_left : this->m_value_right;
  }
};

template <typename T, Curve C, BoundaryCondition BC>
class Extrapolator<T, C, BC, Extrapolation::PERIODIC>
{
private:
  T m_value_left{};
  T m_value_right{};
  T m_period{};

public:
  Extrapolator() = default;

  Extrapolator(Atter<T, C, BC> const &atter, size_t degree)
      : m_value_left{atter.at(degree)}, m_value_right{atter.at(atter.size() - degree - 1)},
        m_period{this->m_value_right - this->m_value_left}
  {
  }

  Extrapolator(Extrapolator const &other) = default;

  Extrapolator(Extrapolator &&other) = default;

  ~Extrapolator() = default;

  Extrapolator &operator=(Extrapolator const &other) = default;

  Extrapolator &operator=(Extrapolator &&other) = default;

  [[nodiscard]] T extrapolate(T value) const
  {
    debugassert(
        value < this->m_value_left or value > this->m_value_right, "Value not outside of the domain"
    );

    T wrapped = std::fmod<T>(value - this->m_value_left, this->m_period);

    if (wrapped < constants::ZERO<T>)
    {
      wrapped += this->m_period;
    }

    return wrapped + this->m_value_left;
  }
};

} // namespace bsplinex::knots

#endif
