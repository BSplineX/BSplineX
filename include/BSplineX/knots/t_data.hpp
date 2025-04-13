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
  T begin{};
  T end{};
  size_t num_elems{0};
  T step_size{};

public:
  Data() { DEBUG_LOG_CALL(); }

  Data(std::vector<T> const &data)
  {
    releaseassert(Data::is_uniform(data), "data is not uniform");

    this->begin     = data.front();
    this->end       = data.back();
    this->num_elems = data.size();
    this->step_size = (this->end - this->begin) / (this->num_elems - 1);
  }

  // Specifying the num-elems means the domain will be [begin, end]
  Data(T begin, T end, size_t num_elems)
  {
    DEBUG_LOG_CALL();
    debugassert(begin < end, "Wrong interval");

    this->begin     = begin;
    this->end       = end;
    this->num_elems = num_elems;
    this->step_size = (end - begin) / (num_elems - 1);
  }

  Data(Data const &other)
      : begin(other.begin), end(other.end), num_elems(other.num_elems), step_size(other.step_size)
  {
    DEBUG_LOG_CALL();
  }

  Data(Data &&other) noexcept
      : begin(other.begin), end(other.end), num_elems(other.num_elems), step_size(other.step_size)
  {
    DEBUG_LOG_CALL();
  }

  ~Data() noexcept { DEBUG_LOG_CALL(); }

  Data &operator=(Data const &other)
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    begin     = other.begin;
    end       = other.end;
    num_elems = other.num_elems;
    step_size = other.step_size;
    return *this;
  }

  Data &operator=(Data &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
    {
      return *this;
    }

    begin     = other.begin;
    end       = other.end;
    num_elems = other.num_elems;
    step_size = other.step_size;
    return *this;
  }

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
    std::generate_n(
        std::back_inserter(tmp), last - first, [this, i = first]() mutable { return this->at(i++); }
    );

    return tmp;
  }

  void pop_tails()
  {
    debugassert(this->num_elems >= 2, "Cannot pop tails from a domain with less than 2 elements");
    this->begin     += this->step_size;
    this->end       -= this->step_size;
    this->num_elems -= 2;
  }

private:
  static bool is_uniform(std::vector<T> const &x)
  {
    if (x.size() < 2)
    {
      return true;
    }

    T const expected_step = std::abs(x.at(1) - x.at(0));

    return std::adjacent_find(
               x.begin(),
               x.end() - 1,
               [expected_step](T a, T b)
               {
                 T const actual_step = std::abs(b - a);
                 T const diff        = std::abs(actual_step - expected_step);
                 T const max_val     = std::max(actual_step, expected_step);

                 return not(diff <= constants::RTOL<T> * max_val or diff <= constants::ATOL<T>);
               }
           ) == x.end() - 1;
  }
};

template <typename T>
class Data<T, Curve::NON_UNIFORM>
{
private:
  std::vector<T> raw_data{};

public:
  Data() { DEBUG_LOG_CALL(); }

  Data(std::vector<T> const &data) : raw_data(data)
  {
    DEBUG_LOG_CALL();
    // NOTE: thank the STL for this wonderful backwards built sort check. Think it as if std::less
    // is <= and std::less_equal is <.
    debugassert(
        std::is_sorted(data.begin(), data.end(), std::less<T>{}),
        "The given data must be sorted respecting the operator <=."
    );
  }

  Data(Data const &other) : raw_data(other.raw_data) { DEBUG_LOG_CALL(); }

  Data(Data &&other) noexcept : raw_data(std::move(other.raw_data)) { DEBUG_LOG_CALL(); }

  ~Data() noexcept { DEBUG_LOG_CALL(); }

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
    DEBUG_LOG_CALL()
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
