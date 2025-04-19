#ifndef BSPLINEX_CONTROL_POINTS_C_PADDER_HPP
#define BSPLINEX_CONTROL_POINTS_C_PADDER_HPP

// Standard includes
#include <vector>

// BSplineX includes
#include "BSplineX/control_points/c_data.hpp"
#include "BSplineX/defines.hpp"
#include "BSplineX/types.hpp"
#include "BSplineX/windows.hpp"

namespace bsplinex::control_points
{

template <typename T, BoundaryCondition BC>
class Padder
{
public:
  Padder() = default;

  Padder(Data<T> & /*data*/, size_t /*degree*/) {}

  Padder(Padder const & /*other*/) = default;

  Padder(Padder && /*other*/) noexcept = default;

  ~Padder() noexcept = default;

  Padder &operator=(Padder const &other) = default;

  Padder &operator=(Padder &&other) noexcept = default;

  [[nodiscard]] T right(size_t /*index*/) const
  {
    releaseassert(
        false,
        "Generic control points padder has zero length, this function is here only for "
        "compatibility reasons."
    );
    return constants::ZERO<T>;
  }

  [[nodiscard]] size_t size() const { return 0; }

  [[nodiscard]] size_t size_right() const { return 0; }
};

template <typename T>
class Padder<T, BoundaryCondition::PERIODIC>
{
private:
  std::vector<T> pad_right{};

public:
  Padder() = default;

  Padder(Data<T> &data, size_t degree) : pad_right{data.slice(0, degree)} {}

  Padder(Padder const &other) = default;

  Padder(Padder &&other) noexcept = default;

  ~Padder() noexcept = default;

  Padder &operator=(Padder const &other) = default;

  Padder &operator=(Padder &&other) noexcept = default;

  [[nodiscard]] T right(size_t index) const
  {
    debugassert(index < this->pad_right.size(), "Out of bounds");
    return this->pad_right.at(index);
  }

  [[nodiscard]] size_t size() const { return this->pad_right.size(); }

  [[nodiscard]] size_t size_right() const { return this->pad_right.size(); }
};

} // namespace bsplinex::control_points

#endif
