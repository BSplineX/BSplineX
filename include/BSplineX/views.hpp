#ifndef BSPLINEX_VIEWS_HPP
#define BSPLINEX_VIEWS_HPP

#include <iterator>

#include "BSplineX/defines.hpp"
#include "BSplineX/windows.hpp"

namespace bsplinex::views
{

/**
 * @brief A view into a contiguous range of elements.
 *
 * This class provides a non-owning view into a range of elements defined by two iterators.
 * It supports basic operations such as element access and size retrieval.
 *
 * @tparam Iter The type of the iterator used to define the range.
 */
template <class Iter>
class ArrayView
{
private:
  Iter _begin;
  Iter _end;

  using difference_type = typename std::iterator_traits<Iter>::difference_type;
  using reference       = typename std::iterator_traits<Iter>::reference;

public:
  /**
   * @brief Default constructor.
   *
   * Constructs an empty ArrayView.
   */
  ArrayView() = default;

  /**
   * @brief Default destructor
   */
  ~ArrayView() = default;

  /**
   * @brief Constructs an ArrayView from a range defined by two iterators.
   *
   * @param begin The iterator pointing to the beginning of the range.
   * @param end The iterator pointing to the end of the range.
   * @throw in debug mode if end is not greater than begin.
   */
  ArrayView(Iter begin, Iter end) : _begin{begin}, _end{end}
  {
    releaseassert(end > begin, "Invalid view range.");
  }

  /**
   * @brief Copy constructor.
   *
   * Constructs an ArrayView by copying another ArrayView.
   *
   * @param other The ArrayView to copy from.
   */
  ArrayView(ArrayView const &other) = default;

  /**
   * @brief Move constructor.
   *
   * Constructs an ArrayView by moving another ArrayView.
   *
   * @param other The ArrayView to move from.
   */
  ArrayView(ArrayView &&other) noexcept = default;

  /**
   * @brief Copy assignment operator.
   *
   * Assigns the contents of another ArrayView to this ArrayView.
   *
   * @param other The ArrayView to copy from.
   * @return A reference to this ArrayView.
   */
  ArrayView &operator=(ArrayView const &other) = default;

  /**
   * @brief Move assignment operator.
   *
   * Assigns the contents of another ArrayView to this ArrayView by moving.
   *
   * @param other The ArrayView to move from.
   * @return A reference to this ArrayView.
   */
  ArrayView &operator=(ArrayView &&other) noexcept = default;

  /**
   * @brief Access an element at a given index with bounds checking.
   *
   * @param index The index of the element to access.
   * @return A reference to the element at the specified index.
   * @throw in debug mode if the index is out of bounds.
   */
  reference at(difference_type index)
  {
    releaseassert(std::distance(this->_begin, this->_end) > index, "Out of bounds.");
    releaseassert(index >= 0, "Negative indices are not supported.");

    return this->operator[](index);
  }

  /**
   * @brief Access an element at a given index without bounds checking.
   *
   * @param index The index of the element to access.
   * @return A reference to the element at the specified index.
   */
  reference operator[](difference_type index) { return *std::next(this->_begin, index); }

  /**
   * @brief Get the first element of the view
   *
   * @return a reference to the first element of the view
   */
  reference front() { return *(this->_begin); }

  /**
   * @brief Get the last element of the view
   *
   * @return a reference to the last element of the view
   */
  reference back() { return *std::prev(this->_end, 1); }

  /**
   * @brief Get the number of elements in the view.
   *
   * @return The number of elements in the view.
   */
  [[nodiscard]] size_t size() const
  {
    return static_cast<size_t>(std::distance(this->_begin, this->_end));
  }

  /**
   * @brief Get the begin iterator
   *
   * @return the begin iterator
   */
  [[nodiscard]] Iter begin() { return this->_begin; }

  /**
   * @brief Get the end iterator
   *
   * @return the end iterator
   */
  [[nodiscard]] Iter end() { return this->_end; }
};

} // namespace bsplinex::views

#endif
