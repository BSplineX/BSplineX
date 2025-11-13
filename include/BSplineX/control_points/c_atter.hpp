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
  Data<T> m_data{};
  Padder<T, BC> m_padder{};

public:
  Atter() = default;

  Atter(Data<T> const &data, size_t degree) : m_data{data}, m_padder{this->m_data, degree} {}

  Atter(Atter const &other) = default;

  Atter(Atter &&other) = default;

  ~Atter() = default;

  Atter &operator=(Atter const &other) = default;

  Atter &operator=(Atter &&other) = default;

  [[nodiscard]] T at(size_t index) const
  {
    debugassert(index < this->size(), "Out of bounds");
    if (index < this->m_data.size())
    {
      return this->m_data.at(index);
    }
    else
    {
      return this->m_padder.right(index - this->m_data.size());
    }
  }

  [[nodiscard]] size_t size() const { return this->m_data.size() + this->m_padder.size(); }

  [[nodiscard]] std::vector<T> get_values() const
  {
    std::vector<T> values;
    values.reserve(m_data.size() + m_padder.size());
    for (size_t i = 0; i < m_data.size(); i++)
    {
      values.push_back(m_data.at(i));
    }
    for (size_t i = 0; i < m_padder.size(); i++)
    {
      values.push_back(m_padder.right(i));
    }
    return values;
  }

  [[nodiscard]] size_t get_derivative_data_size() const
  {
    size_t const data_size = this->m_data.size();
    return this->m_padder.size() == 0 ? data_size - 1 : data_size;
  }
};

} // namespace bsplinex::control_points

#endif
