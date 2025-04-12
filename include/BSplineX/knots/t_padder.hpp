#ifndef BSPLINEX_KNOTS_T_PADDER_HPP
#define BSPLINEX_KNOTS_T_PADDER_HPP

// Standard includes
#include <vector>

// BSplineX includes
#include "BSplineX/defines.hpp"
#include "BSplineX/knots/t_data.hpp"
#include "BSplineX/types.hpp"

namespace bsplinex::knots
{

template <typename T, Curve C, BoundaryCondition BC>
class Padder
{
public:
  [[nodiscard]] virtual T left(size_t index) const  = 0;
  [[nodiscard]] virtual T right(size_t index) const = 0;
  [[nodiscard]] virtual size_t size() const         = 0;
  [[nodiscard]] virtual size_t size_left() const    = 0;
  [[nodiscard]] virtual size_t size_right() const   = 0;
  virtual void pop_tails()                          = 0;
};

template <typename T, Curve C>
class Padder<T, C, BoundaryCondition::OPEN>
{
public:
  Padder() { DEBUG_LOG_CALL(); }

  Padder(Data<T, C> const & /*data*/, size_t /*degree*/) { DEBUG_LOG_CALL(); }

  Padder(Padder const & /*other*/) { DEBUG_LOG_CALL(); }

  Padder(Padder && /*other*/) noexcept { DEBUG_LOG_CALL(); }

  ~Padder() noexcept { DEBUG_LOG_CALL(); }

  Padder &operator=(Padder const &other)
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    return *this;
  }

  Padder &operator=(Padder &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    return *this;
  }

  [[nodiscard]] T left(size_t /*index*/) const
  {
    releaseassert(
        false,
        "OPEN knots padder has zero length, this function is here only for compatibility reasons."
    );
    return constants::ZERO<T>;
  }

  [[nodiscard]] T right(size_t /*index*/) const
  {
    releaseassert(
        false,
        "OPEN knots padder has zero length, this function is here only for compatibility reasons."
    );
    return constants::ZERO<T>;
  }

  [[nodiscard]] size_t size() const { return 0; }

  [[nodiscard]] size_t size_left() const { return 0; }

  [[nodiscard]] size_t size_right() const { return 0; }

  void pop_tails()
  {
    releaseassert(
        false,
        "OPEN knots padder has zero length, this function is here only for compatibility reasons."
    );
  }
};

template <typename T, Curve C>
class Padder<T, C, BoundaryCondition::CLAMPED>
{
private:
  T pad_left{};
  T pad_right{};
  size_t pad_size{0};

public:
  Padder() { DEBUG_LOG_CALL(); }

  Padder(Data<T, C> const &data, size_t degree)
  {
    DEBUG_LOG_CALL();
    this->pad_left  = data.at(0);
    this->pad_right = data.at(data.size() - 1);
    this->pad_size  = degree;
  }

  Padder(Padder const &other)
      : pad_left(other.pad_left), pad_right(other.pad_right), pad_size(other.pad_size)
  {
    DEBUG_LOG_CALL();
  }

  Padder(Padder &&other) noexcept
      : pad_left(other.pad_left), pad_right(other.pad_right), pad_size(other.pad_size)
  {
    DEBUG_LOG_CALL();
  }

  ~Padder() noexcept { DEBUG_LOG_CALL(); }

  Padder &operator=(Padder const &other)
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    pad_left  = other.pad_left;
    pad_right = other.pad_right;
    pad_size  = other.pad_size;
    return *this;
  }

  Padder &operator=(Padder &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    pad_left  = other.pad_left;
    pad_right = other.pad_right;
    pad_size  = other.pad_size;
    return *this;
  }

  [[nodiscard]] T left([[maybe_unused]] size_t index) const
  {
    debugassert(index < this->pad_size, "Out of bounds");
    return this->pad_left;
  }

  [[nodiscard]] T right([[maybe_unused]] size_t index) const
  {
    debugassert(index < this->pad_size, "Out of bounds");
    return this->pad_right;
  }

  [[nodiscard]] size_t size() const { return this->size_left() + this->size_right(); }

  [[nodiscard]] size_t size_left() const { return this->pad_size; }

  [[nodiscard]] size_t size_right() const { return this->pad_size; }

  void pop_tails() { --this->pad_size; }
};

template <typename T, Curve C>
class Padder<T, C, BoundaryCondition::PERIODIC>
{
private:
  std::vector<T> pad_left{};
  std::vector<T> pad_right{};

public:
  Padder() { DEBUG_LOG_CALL(); }

  Padder(Data<T, C> const &data, size_t degree)
  {
    DEBUG_LOG_CALL();
    T period        = data.at(data.size() - 1) - data.at(0);
    this->pad_left  = data.slice(data.size() - degree - 1, data.size() - 1);
    this->pad_right = data.slice(1, degree + 1);
    for (size_t i{0}; i < degree; i++)
    {
      this->pad_left[i]  -= period;
      this->pad_right[i] += period;
    }
  }

  Padder(Padder const &other) : pad_left(other.pad_left), pad_right(other.pad_right)
  {
    DEBUG_LOG_CALL();
  }

  Padder(Padder &&other) noexcept
      : pad_left(std::move(other.pad_left)), pad_right(std::move(other.pad_right))
  {
    DEBUG_LOG_CALL();
  }

  ~Padder() noexcept { DEBUG_LOG_CALL(); }

  Padder &operator=(Padder const &other)
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    pad_left  = other.pad_left;
    pad_right = other.pad_right;
    return *this;
  }

  Padder &operator=(Padder &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    pad_left  = std::move(other.pad_left);
    pad_right = std::move(other.pad_right);
    return *this;
  }

  [[nodiscard]] T left(size_t index) const
  {
    debugassert(index < this->pad_left.size(), "Out of bounds");
    return this->pad_left[index];
  }

  [[nodiscard]] T right(size_t index) const
  {
    debugassert(index < this->pad_right.size(), "Out of bounds");
    return this->pad_right[index];
  }

  [[nodiscard]] size_t size() const { return this->size_left() + this->size_right(); }

  [[nodiscard]] size_t size_left() const { return this->pad_left.size(); }

  [[nodiscard]] size_t size_right() const { return this->pad_right.size(); }

  void pop_tails()
  {
    debugassert(not this->pad_left.empty(), "Cannot pop tails from an empty domain");
    debugassert(not this->pad_right.empty(), "Cannot pop tails from an empty domain");
    this->pad_left.erase(this->pad_left.begin());
    this->pad_right.pop_back();
  }
};

} // namespace bsplinex::knots

#endif
