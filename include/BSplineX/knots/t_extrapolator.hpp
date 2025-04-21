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
  T value_left{};
  T value_right{};

public:
  Extrapolator() = default;

  Extrapolator(Atter<T, C, BC> const &atter, size_t degree)
      : value_left{atter.at(degree)}, value_right{atter.at(atter.size() - degree - 1)}
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
        value < this->value_left or value > this->value_right, "Value not outside of the domain"
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
  Extrapolator() = default;

  Extrapolator(Atter<T, C, BC> const &atter, size_t degree)
      : value_left{atter.at(degree)}, value_right{atter.at(atter.size() - degree - 1)},
        period{this->value_right - this->value_left}
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
        value < this->value_left or value > this->value_right, "Value not outside of the domain"
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
