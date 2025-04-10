#ifndef BSPLINEX_KNOTS_T_EXTRAPOLATOR_HPP
#define BSPLINEX_KNOTS_T_EXTRAPOLATOR_HPP

// Standard includes
#include <cmath>
#include <cstddef>
#include <stdexcept>

// BSplineX includes
#include "BSplineX/defines.hpp"
#include "BSplineX/knots/t_atter.hpp"
#include "BSplineX/types.hpp"

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
  T value_left{};
  T value_right{};

public:
  Extrapolator() { DEBUG_LOG_CALL(); }

  Extrapolator(Atter<T, C, BC> const &atter, size_t degree)
      : value_left{atter.at(degree)}, value_right{atter.at(atter.size() - degree - 1)}
  {
    DEBUG_LOG_CALL();
  }

  Extrapolator(Extrapolator const &other)
      : value_left(other.value_left), value_right(other.value_right)
  {
    DEBUG_LOG_CALL();
  }

  Extrapolator(Extrapolator &&other) noexcept
      : value_left(other.value_left), value_right(other.value_right)
  {
    DEBUG_LOG_CALL();
  }

  ~Extrapolator() { DEBUG_LOG_CALL(); }

  Extrapolator &operator=(Extrapolator const &other)
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    value_left  = other.value_left;
    value_right = other.value_right;
    return *this;
  }

  Extrapolator &operator=(Extrapolator &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    value_left  = other.value_left;
    value_right = other.value_right;
    return *this;
  }

  [[nodiscard]] T extrapolate(T value) const
  {
    debugassert(
        value < this->value_left || value > this->value_right, "Value not outside of the domain"
    );
    return value < this->value_left ? this->value_left : this->value_right;
  }
};

template <typename T, Curve C, BoundaryCondition BC>
class Extrapolator<T, C, BC, Extrapolation::PERIODIC>
{
private:
  T value_left{};
  T value_right{};
  T period{};

public:
  Extrapolator() { DEBUG_LOG_CALL(); }

  Extrapolator(Atter<T, C, BC> const &atter, size_t degree)
      : value_left{atter.at(degree)}, value_right{atter.at(atter.size() - degree - 1)},
        period{this->value_right - this->value_left}
  {
    DEBUG_LOG_CALL();
  }

  Extrapolator(Extrapolator const &other)
      : value_left(other.value_left), value_right(other.value_right), period(other.period)
  {
    DEBUG_LOG_CALL();
  }

  Extrapolator(Extrapolator &&other) noexcept
      : value_left(other.value_left), value_right(other.value_right), period(other.period)
  {
    DEBUG_LOG_CALL();
  }

  ~Extrapolator() { DEBUG_LOG_CALL(); }

  Extrapolator &operator=(Extrapolator const &other)
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    value_left  = other.value_left;
    value_right = other.value_right;
    period      = other.period;
    return *this;
  }

  Extrapolator &operator=(Extrapolator &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    value_left  = other.value_left;
    value_right = other.value_right;
    period      = other.period;
    return *this;
  }

  [[nodiscard]] T extrapolate(T value) const
  {
    debugassert(
        value < this->value_left || value > this->value_right, "Value not outside of the domain"
    );

    T wrapped = std::fmod<T>(value - this->value_left, this->period);

    if (wrapped < constants::ZERO<T>)
    {
      wrapped += this->period;
    }

    return wrapped + this->value_left;
  }
};

} // namespace bsplinex::knots

#endif
