#ifndef T_FINDER_HPP
#define T_FINDER_HPP

// Standard includes
#include <algorithm>
#include <cstddef>
#include <iostream>

// BSplineX includes
#include "BSplineX/defines.hpp"
#include "BSplineX/knots/t_atter.hpp"
#include "BSplineX/types.hpp"

namespace bsplinex::knots
{

template <typename T, Curve C, BoundaryCondition BC, Extrapolation EXT>
class Finder
{
private:
  Atter<T, C, BC> const *atter{nullptr};
  size_t index_left{0};
  size_t index_right{0};

public:
  Finder() = default;

  Finder(Atter<T, C, BC> const &atter, size_t degree)
      : atter{&atter}, index_left{degree}, index_right{this->atter->size() - degree - 1}
  {
  }

  size_t find(T value) const
  {
    // TODO: we have to decide if we want to include the right value of the
    // domain or not
    assertm(
        value >= this->atter->at(this->index_left) && value <= this->atter->at(this->index_right),
        "Value outside of the domain"
    );

    if (!this->atter)
    {
      std::cout << "Helo" << std::endl;
    }
    else
    {
      auto tmp = this->atter->begin();
      std::cout << "Mega helo" << std::endl;
    }

    auto upper = std::upper_bound(
        this->atter->begin() + this->index_left, this->atter->begin() + this->index_right, value
    );

    return upper - this->atter->begin() - 1;
  }
};

template <typename T, BoundaryCondition BC, Extrapolation EXT>
class Finder<T, Curve::UNIFORM, BC, EXT>
{
private:
  T value_left{};
  T value_right{};
  T step_size_inv{};
  size_t degree{};

public:
  Finder() = default;

  Finder(Atter<T, Curve::UNIFORM, BC> const &atter, size_t degree)
      : value_left{atter.at(degree)}, value_right{atter.at(atter.size() - degree - 1)},
        step_size_inv{T(1) / (atter.at(degree + 1) - atter.at(degree))}, degree{degree}
  {
  }

  size_t find(T value) const
  {
    assertm(value >= this->value_left && value <= this->value_right, "Value outside of the domain");

    return static_cast<size_t>((value - this->value_left) * this->step_size_inv) + this->degree;
  }
};

} // namespace bsplinex::knots

#endif