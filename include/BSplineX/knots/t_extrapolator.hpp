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
  virtual size_t extrapolate(T value) const = 0;
};

template <typename T, Curve C, BoundaryCondition BC>
class Extrapolator<T, C, BC, Extrapolation::NONE>
{
public:
  Extrapolator() = default;

  Extrapolator(Atter<T, C, BC> const &, size_t) {}

  T extrapolate(T) const { throw std::runtime_error("Extrapolation explicitly set to NONE"); }
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
      return *this;
    value_left  = other.value_left;
    value_right = other.value_right;
    return *this;
  }

  Extrapolator &operator=(Extrapolator &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
      return *this;
    value_left  = other.value_left;
    value_right = other.value_right;
    return *this;
  }

  T extrapolate(T value) const
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
      return *this;
    value_left  = other.value_left;
    value_right = other.value_right;
    period      = other.period;
    return *this;
  }

  Extrapolator &operator=(Extrapolator &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
      return *this;
    value_left  = other.value_left;
    value_right = other.value_right;
    period      = other.period;
    return *this;
  }

  T extrapolate(T value) const
  {
    debugassert(
        value < this->value_left || value > this->value_right, "Value not outside of the domain"
    );

    // TODO: Figure out how to prevent numerical errors

    if (value < this->value_left)
    {
      value += this->period * (std::floor((this->value_left - value) / this->period) + 1);
    }
    else if (value > this->value_right)
    {
      value -= this->period * (std::floor((value - this->value_right) / this->period) + 1);
    }

    if (value < this->value_left)
    {
      value = this->value_left;
    }
    else if (value > this->value_right)
    {
      value = this->value_right;
    }

    return value;
  }
};

/*

value_right = 1
value = 2
period = 1




*/

} // namespace bsplinex::knots

#endif
