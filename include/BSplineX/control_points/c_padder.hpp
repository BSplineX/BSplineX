#ifndef C_PADDER_HPP
#define C_PADDER_HPP

// Standard includes
#include <stdexcept>
#include <vector>

// BSplineX includes
#include "BSplineX/control_points/c_data.hpp"
#include "BSplineX/defines.hpp"
#include "BSplineX/types.hpp"

namespace bsplinex::control_points
{

template <typename T, BoundaryCondition BC>
class Padder
{
public:
  Padder() { DEBUG_LOG_CALL(); }

  Padder(Data<T> &, size_t) { DEBUG_LOG_CALL(); }

  Padder(Padder const &) { DEBUG_LOG_CALL(); }

  Padder(Padder &&) noexcept { DEBUG_LOG_CALL(); }

  ~Padder() noexcept { DEBUG_LOG_CALL(); }

  Padder &operator=(Padder const &other)
  {
    DEBUG_LOG_CALL();
    if (this == &other)
      return *this;
    return *this;
  }

  Padder &operator=(Padder &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
      return *this;
    return *this;
  }

  T right(size_t) const
  {
    throw std::runtime_error(
        "Generic control points padder has zero length, this function is here only for "
        "compatibility reasons."
    );
  }

  [[nodiscard]] size_t size() const { return 0; }

  [[nodiscard]] size_t size_right() const { return 0; }
};

template <typename T>
class Padder<T, BoundaryCondition::PERIODIC>
{
private:
  std::vector<T> pad_right{};

public:
  Padder() { DEBUG_LOG_CALL(); }

  Padder(Data<T> &data, size_t degree) : pad_right{data.slice(0, degree)} { DEBUG_LOG_CALL(); }

  Padder(Padder const &other) : pad_right(other.pad_right) { DEBUG_LOG_CALL(); }

  Padder(Padder &&other) noexcept : pad_right(std::move(other.pad_right)) { DEBUG_LOG_CALL(); }

  ~Padder() noexcept { DEBUG_LOG_CALL(); }

  Padder &operator=(Padder const &other)
  {
    DEBUG_LOG_CALL();
    if (this == &other)
      return *this;
    pad_right = other.pad_right;
    return *this;
  }

  Padder &operator=(Padder &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
      return *this;
    pad_right = std::move(other.pad_right);
    return *this;
  }

  T right(size_t index) const
  {
    assertm(index < this->pad_right.size(), "Out of bounds");
    return this->pad_right[index];
  }

  [[nodiscard]] size_t size() const { return this->pad_right.size(); }

  [[nodiscard]] size_t size_right() const { return this->pad_right.size(); }
};

} // namespace bsplinex::control_points

#endif
