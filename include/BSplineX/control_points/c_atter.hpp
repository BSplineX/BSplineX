#ifndef BSPLINEX_CONTROL_POINTS_C_ATTER_HPP
#define BSPLINEX_CONTROL_POINTS_C_ATTER_HPP

// BSplineX includes
#include "BSplineX/control_points/c_data.hpp"
#include "BSplineX/control_points/c_padder.hpp"
#include "BSplineX/defines.hpp"
#include "BSplineX/types.hpp"

namespace bsplinex::control_points
{

template <typename T, BoundaryCondition BC>
class Atter
{
private:
  Data<T> data{};
  Padder<T, BC> padder{};

public:
  Atter() { DEBUG_LOG_CALL(); }

  Atter(Data<T> data, size_t degree) : data{data}, padder{this->data, degree} { DEBUG_LOG_CALL(); }

  Atter(Atter const &other) : data(other.data), padder(other.padder) { DEBUG_LOG_CALL(); }

  Atter(Atter &&other) noexcept : data(std::move(other.data)), padder(std::move(other.padder))
  {
    DEBUG_LOG_CALL();
  }

  ~Atter() noexcept { DEBUG_LOG_CALL(); }

  Atter &operator=(Atter const &other)
  {
    DEBUG_LOG_CALL();
    if (this == &other)
      return *this;
    data   = other.data;
    padder = other.padder;
    return *this;
  }

  Atter &operator=(Atter &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
      return *this;
    data   = std::move(other.data);
    padder = std::move(other.padder);
    return *this;
  }

  T at(size_t index) const
  {
    debugassert(index < this->size(), "Out of bounds");
    if (index < this->data.size())
    {
      return this->data.at(index);
    }
    else
    {
      return this->padder.right(index - this->data.size());
    }
  }

  [[nodiscard]] size_t size() const { return this->data.size() + this->padder.size(); }

  std::vector<T> get_values() const
  {
    std::vector<T> values(data.size() + padder.size());
    for (size_t i = 0; i < data.size(); i++)
    {
      values[i] = data.at(i);
    }
    for (size_t i = 0; i < padder.size(); i++)
    {
      values[data.size() + i] = padder.right(i);
    }
    return values;
  }

  [[nodiscard]] size_t get_derivative_data_size() const
  {
    size_t const data_size = this->data.size();
    return (this->padder.size() == 0) ? data_size - 1 : data_size;
  }
};

} // namespace bsplinex::control_points

#endif
