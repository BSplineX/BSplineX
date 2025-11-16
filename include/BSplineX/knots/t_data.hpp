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
#include "BSplineX/windows.hpp"

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
  T m_begin{};
  T m_end{};
  size_t m_num_elems{0};
  T m_step_size{};

public:
  Data() = default;

  explicit Data(std::vector<T> const &data)
  {
    releaseassert(Data::is_uniform(data), "Data must be uniform with step > 0");

    this->m_begin     = data.front();
    this->m_end       = data.back();
    this->m_num_elems = data.size();
    this->m_step_size = (this->m_end - this->m_begin) / (this->m_num_elems - 1);
  }

  // Specifying the num-elems means the domain will be [begin, end]
  Data(T begin, T end, size_t num_elems)
  {
    debugassert(begin < end, "Wrong interval");

    this->m_begin     = begin;
    this->m_end       = end;
    this->m_num_elems = num_elems;
    this->m_step_size = (end - begin) / (num_elems - 1);
  }

  Data(Data const &other) = default;

  Data(Data &&other) = default;

  ~Data() = default;

  Data &operator=(Data const &other) = default;

  Data &operator=(Data &&other) = default;

  [[nodiscard]] T at(size_t index) const
  {
    debugassert(index < this->m_num_elems, "Out of bounds");
    return std::fma<T>(static_cast<T>(index), this->m_step_size, this->m_begin);
  }

  [[nodiscard]] size_t size() const { return this->m_num_elems; }

  [[nodiscard]] std::vector<T> slice(size_t first, size_t last) const
  {
    debugassert(first <= last, "Invalid range");
    debugassert(last <= this->m_num_elems, "Out of bounds");

    std::vector<T> tmp{};
    tmp.reserve(last - first);
    std::generate_n(
        std::back_inserter(tmp), last - first, [this, i = first]() mutable { return this->at(i++); }
    );

    return tmp;
  }

  void pop_tails()
  {
    debugassert(this->m_num_elems >= 2, "Cannot pop tails from a domain with less than 2 elements");
    this->m_begin     += this->m_step_size;
    this->m_end       -= this->m_step_size;
    this->m_num_elems -= 2;
  }

private:
  static bool is_uniform(std::vector<T> const &x)
  {
    if (x.size() < 2)
    {
      return true;
    }

    T const expected_step = x.at(1) - x.at(0);
    if (expected_step <= 0)
    {
      return false;
    }

    return std::adjacent_find(
               x.begin(),
               x.end(),
               [expected_step](T a, T b)
               {
                 T const actual_step = b - a;
                 T const diff        = std::abs(actual_step - expected_step);
                 T const max_val     = std::max(actual_step, expected_step);

                 return not(diff <= constants::RTOL<T> * max_val or diff <= constants::ATOL<T>);
               }
           ) == x.end();
  }
};

template <typename T>
class Data<T, Curve::NON_UNIFORM>
{
private:
  std::vector<T> m_raw_data{};

public:
  Data() = default;

  explicit Data(std::vector<T> const &data) : m_raw_data(data)
  {
    // NOTE: thank the STL for this wonderful backwards built sort check. Think it as if std::less
    // is <= and std::less_equal is <.
    debugassert(
        std::is_sorted(data.begin(), data.end(), std::less<T>{}),
        "The given data must be sorted respecting the operator <=."
    );
  }

  Data(Data const &other) = default;

  Data(Data &&other) = default;

  ~Data() = default;

  Data &operator=(Data const &other) = default;

  Data &operator=(Data &&other) = default;

  [[nodiscard]] T at(size_t index) const
  {
    debugassert(index < this->m_raw_data.size(), "Out of bounds");
    return this->m_raw_data[index];
  }

  [[nodiscard]] size_t size() const { return this->m_raw_data.size(); }

  [[nodiscard]] std::vector<T> slice(size_t first, size_t last) const
  {
    debugassert(first <= last, "Invalid range");
    debugassert(last <= this->m_raw_data.size(), "Out of bounds");

    using difference_type = typename std::vector<T>::iterator::difference_type;

    return std::vector<T>{
        std::next(this->m_raw_data.begin(), static_cast<difference_type>(first)),
        std::next(this->m_raw_data.begin(), static_cast<difference_type>(last))
    };
  }

  void pop_tails()
  {
    debugassert(
        this->m_raw_data.size() >= 2, "Cannot pop tails from a domain with less than 2 elements"
    );
    this->m_raw_data.pop_back();
    this->m_raw_data.erase(this->m_raw_data.begin());
  }
};

} // namespace bsplinex::knots

#endif
