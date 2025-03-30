#ifndef BSPLINEX_CONTROL_POINTS_C_DATA_HPP
#define BSPLINEX_CONTROL_POINTS_C_DATA_HPP

// Standard includes
#include <cstddef>
#include <vector>

// BSplineX includes
#include "BSplineX/defines.hpp"

namespace bsplinex::control_points
{

template <typename T>
class Data
{
private:
  std::vector<T> raw_data{};

public:
  Data() { DEBUG_LOG_CALL(); }

  Data(std::vector<T> const &data) : raw_data{data} { DEBUG_LOG_CALL(); }

  Data(Data const &other) : raw_data(other.raw_data) { DEBUG_LOG_CALL(); }

  Data(Data &&other) noexcept : raw_data(std::move(other.raw_data)) { DEBUG_LOG_CALL(); }

  ~Data() { DEBUG_LOG_CALL(); }

  Data &operator=(Data const &other)
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    raw_data = other.raw_data;
    return *this;
  }

  Data &operator=(Data &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    raw_data = std::move(other.raw_data);
    return *this;
  }

  [[nodiscard]] T at(size_t index) const
  {
    debugassert(index < this->raw_data.size(), "Out of bounds");
    return this->raw_data[index];
  }

  [[nodiscard]] size_t size() const { return this->raw_data.size(); }

  std::vector<T> slice(size_t first, size_t last)
  {
    debugassert(first <= last, "Invalid range");
    debugassert(last <= this->raw_data.size(), "Out of bounds");

    return std::vector<T>{this->raw_data.begin() + first, this->raw_data.begin() + last};
  }
};

} // namespace bsplinex::control_points

#endif
