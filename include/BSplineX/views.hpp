#ifndef BSPLINEX_VIEWS_HPP
#define BSPLINEX_VIEWS_HPP

#include <iterator>

#include "BSplineX/defines.hpp"

namespace bsplinex::views
{

template <class Iter>
class ArrayView
{
private:
  Iter begin;
  Iter end;

  using difference_type = typename std::iterator_traits<Iter>::difference_type;
  using reference       = typename std::iterator_traits<Iter>::reference;

public:
  ArrayView() {}

  ArrayView(Iter begin, Iter end) : begin{begin}, end{end}
  {
    debugassert(end > begin, "Invalid view range.");
  }

  ArrayView(ArrayView const &other) : begin{other.begin}, end{other.end} {}

  ArrayView(ArrayView &&other) : begin{std::move(other.begin)}, end{std::move(other.end)} {}

  ArrayView &operator=(ArrayView const &other)
  {
    this->begin = other.begin;
    this->end   = other.end;
    return *this;
  }

  ArrayView &operator=(ArrayView &&other)
  {
    this->begin = std::move(other.begin);
    this->end   = std::move(other.end);
    return *this;
  }

  reference at(difference_type index)
  {
    debugassert(std::distance(this->begin, this->end) > index, "Out of bounds.");
    debugassert(index >= 0, "Negative indices are not supported.");

    return this->operator[](index);
  }

  reference const at(difference_type index) const
  {
    debugassert(std::distance(this->begin, this->end) > index, "Out of bounds.");
    debugassert(index >= 0, "Negative indices are not supported.");

    return this->operator[](index);
  }

  reference operator[](difference_type index) { return *std::next(this->begin, index); }

  reference const operator[](difference_type index) const { return *std::next(this->begin, index); }

  size_t size() const { return static_cast<size_t>(std::distance(this->begin, this->end)); }
};

} // namespace bsplinex::views

#endif
