#ifndef BSPLINEX_CONTROL_POINTS_C_ATTER_HPP
#define BSPLINEX_CONTROL_POINTS_C_ATTER_HPP

// BSplineX includes
#include "BSplineX/control_points/c_data.hpp"
#include "BSplineX/control_points/c_padder.hpp"
#include "BSplineX/defines.hpp"
#include "BSplineX/types.hpp"
#include "BSplineX/windows.hpp"

namespace bsplinex::control_points
{

template <typename T, BoundaryCondition BC>
class Atter
{
private:
  Data<T> data{};
  Padder<T, BC> padder{};

public:
  Atter() = default;

  Atter(Data<T> data, size_t degree) : data{std::move(data)}, padder{this->data, degree} {}

  Atter(Atter const &other) = default;

  Atter(Atter &&other) noexcept = default;

  ~Atter() noexcept = default;

  Atter &operator=(Atter const &other) = default;

  Atter &operator=(Atter &&other) noexcept = default;

  [[nodiscard]] T at(size_t index) const
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

  [[nodiscard]] std::vector<T> get_values() const
  {
    std::vector<T> values;
    values.reserve(data.size() + padder.size());
    for (size_t i = 0; i < data.size(); i++)
    {
      values.push_back(data.at(i));
    }
    for (size_t i = 0; i < padder.size(); i++)
    {
      values.push_back(padder.right(i));
    }
    return values;
  }

  [[nodiscard]] size_t get_derivative_data_size() const
  {
    size_t const data_size = this->data.size();
    return this->padder.size() == 0 ? data_size - 1 : data_size;
  }
};

} // namespace bsplinex::control_points

#endif
