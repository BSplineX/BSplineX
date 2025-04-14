#ifndef BSPLINEX_KNOTS_T_DATA_HPP
#define BSPLINEX_KNOTS_T_DATA_HPP

// Standard includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

// BSplineX includes
#include "BSplineX/defines.hpp"
#include "BSplineX/types.hpp"

namespace bsplinex::knots
{

template <typename T, Curve C>
class Data
{
public:
  [[nodiscard]] virtual T at(size_t index) const                              = 0;
  [[nodiscard]] virtual size_t size() const                                   = 0;
  [[nodiscard]] virtual std::vector<T> slice(size_t first, size_t last) const = 0;
  virtual void pop_tails()                                                    = 0;
};

template <typename T>
class Data<T, Curve::UNIFORM>
{
private:
  T begin{};
  T end{};
  size_t num_elems{0};
  T step_size{};

public:
  Data() = default;

  // Specifying the step-size means the domain will be [begin, end[
  Data(T begin, T end, T step)
  {
    debugassert(step > 0, "Negative step-size");
    debugassert(begin < end, "Wrong interval");

    this->begin     = begin;
    this->step_size = step;
    this->num_elems = (end - begin) / step + 1;
    this->end       = begin + (this->num_elems - 1) * step;
  }

  // Specifying the num-elems means the domain will be [begin, end]
  Data(T begin, T end, size_t num_elems)
  {
    debugassert(begin < end, "Wrong interval");

    this->begin     = begin;
    this->end       = end;
    this->num_elems = num_elems;
    this->step_size = (end - begin) / (num_elems - 1);
  }

  Data(Data const &other) = default;

  Data(Data &&other) noexcept = default;

  ~Data() noexcept = default;

  Data &operator=(Data const &other) = default;

  Data &operator=(Data &&other) noexcept = default;

  [[nodiscard]] T at(size_t index) const
  {
    debugassert(index < this->num_elems, "Out of bounds");
    return std::fma<T>(static_cast<T>(index), this->step_size, this->begin);
  }

  [[nodiscard]] size_t size() const { return this->num_elems; }

  [[nodiscard]] std::vector<T> slice(size_t first, size_t last) const
  {
    debugassert(first <= last, "Invalid range");
    debugassert(last <= this->num_elems, "Out of bounds");

    std::vector<T> tmp{};
    tmp.reserve(last - first);

    for (size_t i{first}; i < last; i++)
    {
      tmp.push_back(this->at(i));
    }

    return tmp;
  }

  void pop_tails()
  {
    debugassert(this->num_elems >= 2, "Cannot pop tails from a domain with less than 2 elements");
    this->begin     += this->step_size;
    this->end       -= this->step_size;
    this->num_elems -= 2;
  }
};

template <typename T>
class Data<T, Curve::NON_UNIFORM>
{
private:
  std::vector<T> raw_data{};

public:
  Data() = default;

  Data(std::vector<T> const &data) : raw_data(data)
  {
    // NOTE: thank the STL for this wonderful backwards built sort check. Think it as if std::less
    // is <= and std::less_equal is <.
    debugassert(
        std::is_sorted(data.begin(), data.end(), std::less<T>{}),
        "The given data must be sorted respecting the operator <=."
    );
  }

  Data(Data const &other) = default;

  Data(Data &&other) noexcept = default;

  ~Data() noexcept = default;

  Data &operator=(Data const &other) = default;

  Data &operator=(Data &&other) noexcept = default;

  [[nodiscard]] T at(size_t index) const
  {
    debugassert(index < this->raw_data.size(), "Out of bounds");
    return this->raw_data[index];
  }

  [[nodiscard]] size_t size() const { return this->raw_data.size(); }

  [[nodiscard]] std::vector<T> slice(size_t first, size_t last) const
  {
    debugassert(first <= last, "Invalid range");
    debugassert(last <= this->raw_data.size(), "Out of bounds");

    return std::vector<T>{this->raw_data.begin() + first, this->raw_data.begin() + last};
  }

  void pop_tails()
  {
    debugassert(
        this->raw_data.size() >= 2, "Cannot pop tails from a domain with less than 2 elements"
    );
    this->raw_data.pop_back();
    this->raw_data.erase(this->raw_data.begin());
  }
};

} // namespace bsplinex::knots

#endif
