#ifndef BSPLINEX_KNOTS_T_ATTER_HPP
#define BSPLINEX_KNOTS_T_ATTER_HPP

// BSplineX includes
#include "BSplineX/defines.hpp"
#include "BSplineX/knots/t_data.hpp"
#include "BSplineX/knots/t_padder.hpp"
#include "BSplineX/types.hpp"

namespace bsplinex::knots
{

template <typename T, Curve C, BoundaryCondition BC>
class Atter
{
private:
  Data<T, C> data{};
  Padder<T, C, BC> padder{};

public:
  Atter() { DEBUG_LOG_CALL(); }

  Atter(Data<T, C> const &data, size_t degree) : data{data}, padder{this->data, degree}
  {
    DEBUG_LOG_CALL();
  }

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
    if (index < this->padder.size_left())
    {
      return this->padder.left(index);
    }
    else if (index > this->data.size() - 1 + this->padder.size_left())
    {
      return this->padder.right(index - this->data.size() - this->padder.size_left());
    }
    else
    {
      return this->data.at(index - this->padder.size_left());
    }
  }

  [[nodiscard]] size_t size() const { return this->data.size() + this->padder.size(); }

  Atter &pop_tails()
  {
    if (this->padder.size() > 0)
    {
      this->padder.pop_tails();
    }
    else
    {
      this->data.pop_tails();
    }

    return *this;
  }

  class iterator
  {
  private:
    Atter const *atter{nullptr};
    size_t index{0};

  public:
    // iterator traits
    using difference_type   = int;
    using value_type        = T;
    using pointer           = const T *;
    using reference         = const T &;
    using iterator_category = std::random_access_iterator_tag;

    iterator(Atter<T, C, BC> const *atter, size_t index) : atter{atter}, index{index} {}

    iterator(iterator const &b) = default;

    iterator &operator++()
    {
      ++(this->index);
      return *this;
    }

    iterator operator++(int)
    {
      iterator retval = *this;
      ++(*this);
      return retval;
    }

    iterator &operator--()
    {
      --(this->index);
      return *this;
    }

    iterator operator--(int)
    {
      iterator retval = *this;
      --(*this);
      return retval;
    }

    iterator &operator+=(int n)
    {
      this->index += n;
      return *this;
    }

    iterator operator+(int n) const
    {
      iterator retval  = *this;
      retval          += n;
      return retval;
    }

    iterator &operator-=(int n)
    {
      this->index -= n;
      return *this;
    }

    iterator operator-(int n) const
    {
      iterator retval  = *this;
      retval          -= n;
      return retval;
    }

    difference_type operator-(iterator const &b) const { return this->index - b.index; }

    bool operator==(iterator const &other) const { return this->index == other.index; }

    iterator &operator=(iterator const &b)
    {
      if (this == &b)
      {
        return *this;
      }

      this->atter = b.atter;
      this->index = b.index;
      return *this;
    };

    bool operator!=(iterator const &other) const { return !(*this == other); }

    value_type operator*() const { return this->atter->at(this->index); }

    value_type operator[](int n) const { return *(*this + n); }

    bool operator<(iterator const &b) const { return this->index < b.index; }

    bool operator>(iterator const &b) const { return this->index > b.index; }

    bool operator<=(iterator const &b) const { return !(*this > b); }

    bool operator>=(iterator const &b) const { return !(*this < b); }
  };

  iterator begin() const { return {this, 0}; }

  iterator end() const { return {this, this->size()}; }
};

} // namespace bsplinex::knots

#endif
