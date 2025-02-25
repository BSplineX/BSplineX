#ifndef BSPLINEX_KNOTS_T_FINDER_HPP
#define BSPLINEX_KNOTS_T_FINDER_HPP

// Standard includes
#include <algorithm>
#include <cstddef>

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
  Finder() { DEBUG_LOG_CALL(); }

  Finder(Atter<T, C, BC> const &atter, size_t degree)
      : atter{&atter}, index_left{degree}, index_right{this->atter->size() - degree - 1}
  {
    DEBUG_LOG_CALL();
  }

  Finder(Finder const &other) = delete;

  Finder(Finder &&other) = delete;

  ~Finder() noexcept { DEBUG_LOG_CALL(); }

  Finder &operator=(Finder const &other) = delete;

  Finder &operator=(Finder &&other) = delete;

  size_t find(T value) const
  {
    debugassert(
        value >= this->atter->at(this->index_left) && value <= this->atter->at(this->index_right),
        "Value outside of the domain"
    );

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
  size_t max_index{};

public:
  Finder() { DEBUG_LOG_CALL(); }

  Finder(Atter<T, Curve::UNIFORM, BC> const &atter, size_t degree)
      : value_left{atter.at(degree)}, value_right{atter.at(atter.size() - 1 - degree)},
        step_size_inv{static_cast<T>(1) / (atter.at(degree + 1) - atter.at(degree))},
        degree{degree}, max_index{atter.size() - 1 - degree - 1}
  {
    DEBUG_LOG_CALL();
  }

  Finder(Finder const &other) = delete;

  Finder(Finder &&other) = delete;

  Finder &operator=(Finder const &other) = delete;

  Finder &operator=(Finder &&other) = delete;

  size_t find(T value) const
  {
    debugassert(
        value >= this->value_left && value <= this->value_right, "Value outside of the domain"
    );

    size_t index =
        static_cast<size_t>((value - this->value_left) * this->step_size_inv) + this->degree;
    return std::min(index, max_index);
  }
};

// t_i = step_size * i + start_t
//
/*

(begin + (atter.size - 1 - degree) * step_size - begin - degree * step_size) / step_size + degree
((atter.size - 1 - degree) - degree) + degree
atter.size - 1 - degree

knots - 1 - degree = ctrl

*/

} // namespace bsplinex::knots

#endif
