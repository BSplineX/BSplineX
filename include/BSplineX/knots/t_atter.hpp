#ifndef BSPLINEX_KNOTS_T_ATTER_HPP
#define BSPLINEX_KNOTS_T_ATTER_HPP

// Standard includes
#include <cstddef>

// BSplineX includes
#include "BSplineX/defines.hpp"
#include "BSplineX/knots/t_data.hpp"
#include "BSplineX/knots/t_padder.hpp"
#include "BSplineX/types.hpp"
#include "BSplineX/windows.hpp"

namespace bsplinex::knots
{

template <typename T, Curve C, BoundaryCondition BC>
class Atter
{
private:
  Data<T, C> m_data{};
  Padder<T, C, BC> m_padder{};

public:
  Atter() = default;

  Atter(Data<T, C> const &data, size_t degree) : m_data{data}, m_padder{this->m_data, degree} {}

  Atter(Atter const &other) = default;

  Atter(Atter &&other) = default;

  ~Atter() = default;

  Atter &operator=(Atter const &other) = default;

  Atter &operator=(Atter &&other) = default;

  [[nodiscard]] T at(size_t index) const
  {
    debugassert(index < this->size(), "Out of bounds");
    if (index < this->m_padder.size_left())
    {
      return this->m_padder.left(index);
    }
    else if (index > this->m_data.size() - 1 + this->m_padder.size_left())
    {
      return this->m_padder.right(index - this->m_data.size() - this->m_padder.size_left());
    }
    else
    {
      return this->m_data.at(index - this->m_padder.size_left());
    }
  }

  [[nodiscard]] size_t size() const { return this->m_data.size() + this->m_padder.size(); }

  Atter &pop_tails()
  {
    if (this->m_padder.size() > 0)
    {
      this->m_padder.pop_tails();
    }
    else
    {
      this->m_data.pop_tails();
    }

    return *this;
  }

  class iterator
  {
  private:
    Atter const *m_atter{nullptr};
    size_t m_index{0};

  public:
    // iterator traits
    using difference_type   = std::ptrdiff_t;
    using value_type        = T;
    using pointer           = T const *;
    using reference         = T const &;
    using iterator_category = std::random_access_iterator_tag;

    iterator(Atter<T, C, BC> const *atter, size_t index) : m_atter{atter}, m_index{index} {}

    ~iterator() = default;

    iterator(iterator const &b) = default;

    iterator(iterator &&b) = default;

    iterator &operator=(iterator const &b) = default;

    iterator &operator=(iterator &&b) = default;

    iterator &operator++()
    {
      ++(this->m_index);
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
      --(this->m_index);
      return *this;
    }

    iterator operator--(int)
    {
      iterator retval = *this;
      --(*this);
      return retval;
    }

    iterator &operator+=(difference_type n)
    {
      this->m_index += n;
      return *this;
    }

    iterator operator+(difference_type n) const
    {
      iterator retval  = *this;
      retval          += n;
      return retval;
    }

    iterator &operator-=(difference_type n)
    {
      this->m_index -= n;
      return *this;
    }

    iterator operator-(difference_type n) const
    {
      iterator retval  = *this;
      retval          -= n;
      return retval;
    }

    difference_type operator-(iterator const &b) const
    {
      return static_cast<difference_type>(this->m_index - b.m_index);
    }

    bool operator==(iterator const &other) const { return this->m_index == other.m_index; }

    bool operator!=(iterator const &other) const { return !(*this == other); }

    value_type operator*() const { return this->m_atter->at(this->m_index); }

    value_type operator[](difference_type n) const { return *(*this + n); }

    bool operator<(iterator const &b) const { return this->m_index < b.m_index; }

    bool operator>(iterator const &b) const { return this->m_index > b.m_index; }

    bool operator<=(iterator const &b) const { return !(*this > b); }

    bool operator>=(iterator const &b) const { return !(*this < b); }
  };

  [[nodiscard]] iterator begin() const { return {this, 0}; }

  [[nodiscard]] iterator end() const { return {this, this->size()}; }
};

} // namespace bsplinex::knots

#endif
