#ifndef BSPLINEX_BSPLINE_BSPLINE_HPP
#define BSPLINEX_BSPLINE_BSPLINE_HPP

// Standard includes
#include <algorithm>
#include <limits>
#include <sstream>
#include <vector>

// BSplineX includes
#include "BSplineX/bspline/bspline_lsq.hpp"
#include "BSplineX/control_points/control_points.hpp"
#include "BSplineX/defines.hpp"
#include "BSplineX/knots/knots.hpp"
#include "BSplineX/types.hpp"

namespace bsplinex::bspline
{

/**
 * @brief A BSpline class.
 *
 * @tparam T The type of the control points.
 * @tparam C The curve type.
 * @tparam BC The boundary condition.
 * @tparam EXT The extrapolation type.
 */
template <typename T, Curve C, BoundaryCondition BC, Extrapolation EXT>
class BSpline
{
private:
  using vec_iter = typename std::vector<T>::const_iterator;
  using vec_view = views::ArrayView<vec_iter>;

  knots::Knots<T, C, BC, EXT> knots{};
  control_points::ControlPoints<T, BC> control_points{};
  size_t degree{};
  std::vector<T> support{};

public:
  /**
   * @brief Default constructor.
   *
   * Constructs an empty BSpline.
   */
  BSpline() { DEBUG_LOG_CALL(); }

  /**
   * @brief Constructs a BSpline from a knot vector and control points.
   *
   * Note that depending on the boundary condition, the knot vector and control
   * points will be padded to respect the boundary condition.
   *
   * @param knots_data The knot vector data.
   * @param control_points_data The control points data.
   * @param degree The degree of the BSpline.
   */
  BSpline(
      knots::Data<T, C> const &knots_data,
      control_points::Data<T> const &control_points_data,
      size_t degree
  )
      : knots{knots_data, degree}, control_points{control_points_data, degree}, degree{degree}
  {
    DEBUG_LOG_CALL();
    this->check_sizes();
    this->support.resize(this->degree + 1);
  }

  /**
   * @brief Copy constructor.
   *
   * Constructs a BSpline by copying another BSpline.
   *
   * @param other The BSpline to copy from.
   */
  BSpline(BSpline const &other)
      : knots(other.knots), control_points(other.control_points), degree(other.degree),
        support(other.support)
  {
    DEBUG_LOG_CALL();
  }

  /**
   * @brief Move constructor.
   *
   * Constructs a BSpline by moving another BSpline.
   *
   * @param other The BSpline to move from.
   */
  BSpline(BSpline &&other) noexcept
      : knots(std::move(other.knots)), control_points(std::move(other.control_points)),
        degree(other.degree), support(std::move(other.support))
  {
    DEBUG_LOG_CALL();
  }

  /**
   * @brief Destructor.
   */
  ~BSpline() noexcept { DEBUG_LOG_CALL(); }

  /**
   * @brief Copy assignment operator.
   *
   * Assigns the contents of another BSpline to this BSpline.
   *
   * @param other The BSpline to copy from.
   * @return A reference to this BSpline.
   */
  BSpline &operator=(BSpline const &other)
  {
    DEBUG_LOG_CALL();
    if (this == &other)
      return *this;
    knots          = other.knots;
    control_points = other.control_points;
    degree         = other.degree;
    support        = other.support;
    return *this;
  }

  /**
   * @brief Move assignment operator.
   *
   * Assigns the contents of another BSpline to this BSpline by moving.
   *
   * @param other The BSpline to move from.
   * @return A reference to this BSpline.
   */
  BSpline &operator=(BSpline &&other) noexcept
  {
    DEBUG_LOG_CALL();
    if (this == &other)
      return *this;
    knots          = std::move(other.knots);
    control_points = std::move(other.control_points);
    degree         = other.degree;
    support        = std::move(other.support);
    return *this;
  }

  /**
   * @brief Evaluates the BSpline at a given value.
   *
   * @param value The value to evaluate the BSpline at.
   * @return The value of the BSpline at the given value.
   */
  T evaluate(T value)
  {
    auto index_value_pair = this->knots.find(value);
    return this->deboor(index_value_pair.first, index_value_pair.second);
  }

  /**
   * @brief Computes the basis functions at a given value.
   *
   * @param value The value to compute the basis functions at.
   * @return The basis functions at the given value.
   */
  std::vector<T> basis(T value)
  {
    std::vector<T> basis_functions(this->degree + 1, (T)0);

    size_t index = this->nnz_basis(value, basis_functions.begin(), basis_functions.end());

    basis_functions.insert(basis_functions.begin(), index, (T)0);
    basis_functions.insert(
        basis_functions.end(), this->control_points.size() - index - this->degree - 1, (T)0
    );

    return basis_functions;
  }

  /**
   * @brief Fits the BSpline to some data.
   *
   * @param x The x values of the data.
   * @param y The y values of the data.
   */
  void fit(std::vector<T> const &x, std::vector<T> const &y)
  {
    releaseassert(x.size() == y.size(), "x and y must have the same size");

    this->control_points = std::move(lsq::lsq<T, vec_iter, BC>(
        degree,
        knots.size(),
        [this](T value, std::vector<T> &vec) -> size_t
        { return this->nnz_basis(value, vec.begin(), vec.end()); },
        vec_view{x.begin(), x.end()},
        vec_view{y.begin(), y.end()},
        {}
    ));
  }

  /**
   * @brief Interpolates the BSpline to some data.
   *
   * Note that if the boundary condition is:
   * - OPEN: the first and last `degree` x and y values are not interpolated. If
   * you want to use them, either use another boundary condition or manually pad
   * the data.
   * - CLAMPED: all x and y values are interpolated.
   * - PERIODIC: the last x and y values are not interpolated. This seems to be
   * an issue, but if your data is truly periodic, you can simply repeat the
   * first x and y values at the end of the data and everything will work as
   * expected.
   *
   * @param x The x values of the data. They must be sorted in ascending order.
   * @param y The y values of the data.
   * @param additional_conditions Additional conditions to impose on the BSpline.
   * They must lie inside the knots interval. The number of additional conditions
   * must be equal to `degree - 1` except for periodic BSplines where no
   * additional conditions can be provided.
   */
  void interpolate(
      std::vector<T> const &x,
      std::vector<T> const &y,
      std::vector<lsq::Condition<T>> const &additional_conditions
  )
  {
    releaseassert(x.size() == y.size(), "x and y must have the same size");

    if constexpr (BoundaryCondition::PERIODIC == BC)
    {
      releaseassert(
          additional_conditions.empty(),
          "For PERIODIC BSplines there must be exactly 0 additional conditions."
      );
    }
    else
    {
      releaseassert(
          additional_conditions.size() == degree - 1,
          "There must be exactly degree - 1 additional conditions."
      );
    }

    knots::Knots<T, C, BC, EXT> new_knots;
    if constexpr (Curve::UNIFORM == C)
    {
      T step = x.at(1) - x.at(0);
      for (size_t i{0}; i < x.size() - 1; i++)
      {
        releaseassert(
            std::abs(x.at(i + 1) - x.at(i) - step) <= std::numeric_limits<T>::epsilon(),
            "x is not uniform."
        );
      }
      new_knots = std::move(knots::Knots<T, C, BC, EXT>{{x.front(), x.back(), x.size()}, degree});
    }
    else
    {
      new_knots = std::move(knots::Knots<T, C, BC, EXT>{{x}, degree});
    }

    auto const domain = new_knots.domain();
    releaseassert(
        std::all_of(
            additional_conditions.begin(),
            additional_conditions.end(),
            [&domain](auto const &elem)
            { return elem.x_value >= domain.first and elem.x_value <= domain.second; }
        ),
        "Additional conditions must lie inside the knots interval."
    );

    this->knots = std::move(new_knots);

    vec_view x_view, y_view;
    if constexpr (BoundaryCondition::OPEN == BC)
    {
      x_view = vec_view{std::next(x.begin(), degree), std::prev(x.end(), degree)};
      y_view = vec_view{std::next(y.begin(), degree), std::prev(y.end(), degree)};
    }
    else if constexpr (BoundaryCondition::CLAMPED == BC)
    {
      x_view = vec_view{x.begin(), x.end()};
      y_view = vec_view{y.begin(), y.end()};
    }
    else if constexpr (BoundaryCondition::PERIODIC == BC)
    {
      x_view = vec_view{x.begin(), std::prev(x.end(), 1)};
      y_view = vec_view{y.begin(), std::prev(y.end(), 1)};
    }
    else
    {
      releaseassert(false, "Unkown boundary condition, you should never get here!");
    }

    this->control_points = std::move(lsq::lsq<T, vec_iter, BC>(
        degree,
        this->knots.size(),
        [this](T value, std::vector<T> &vec) -> size_t
        { return this->nnz_basis(value, vec.begin(), vec.end()); },
        x_view,
        y_view,
        additional_conditions
    ));
  }

  std::vector<T> get_control_points()
  {
    std::vector<T> ctrl_pts{};
    ctrl_pts.reserve(this->control_points.size());

    for (size_t i{0}; i < this->control_points.size(); i++)
    {
      ctrl_pts.push_back(this->control_points.at(i));
    }

    return ctrl_pts;
  }

private:
  void check_sizes()
  {
    if (this->control_points.size() == this->knots.size() - this->degree - 1)
    {
      return;
    }

    std::stringstream ss{};
    ss << "Found control_points.size() != knots.size() - degree - 1 ("
       << this->control_points.size() << " != " << this->knots.size() - this->degree - 1 << "). ";

    // clang-format off

    if constexpr (BC == BoundaryCondition::OPEN)
    {
      ss << "With BoundaryCondition::OPEN no padding is added, therefore you need to respect: control_points_data.size() = knots_data.size() - degree - 1";
    }
    else if constexpr (BC == BoundaryCondition::CLAMPED)
    {
      ss << "With BoundaryCondition::CLAMPED padding is added to the knots, therefore you need to respect: control_points_data.size() = knots_data.size() + degree - 1";
    }
    else if constexpr (BC == BoundaryCondition::PERIODIC)
    {
      ss << "With BoundaryCondition::PERIODIC padding is added to the knots and control points, therefore you need to respect: control_points_data.size() = knots_data.size() - 1";
    }
    else {
     ss << "Unknown BoundaryCondition, you should not have arrived here ever!";
    }

    // clang-format on

    releaseassert(false, ss.str());
  }

  template <typename It>
  size_t nnz_basis(T value, [[maybe_unused]] It begin, It end)
  {
    debugassert(
        (end - begin) == (long long)(this->degree + 1),
        "Unexpected number of basis asked, exactly degree + 1 basis can be asked"
    );

    debugassert(
        std::all_of(begin, end, [](T i) { return (T)0 == i; }),
        "Initial basis must be initialised to zero"
    );

    auto [index, val] = this->knots.find(value);

    *(end - 1) = 1.0;
    for (size_t d{1}; d <= this->degree; d++)
    {
      *(end - 1 - d) = (this->knots.at(index + 1) - val) /
                       (this->knots.at(index + 1) - this->knots.at(index - d + 1)) *
                       *(end - 1 - d + 1);
      for (size_t i{index - d + 1}; i < index; i++)
      {
        *(end - 1 - index + i) =
            (val - this->knots.at(i)) / (this->knots.at(i + d) - this->knots.at(i)) *
                *(end - 1 - index + i) +
            (this->knots.at(i + d + 1) - val) /
                (this->knots.at(i + d + 1) - this->knots.at(i + 1)) * *(end - 1 - index + i + 1);
      }
      *(end - 1) = (val - this->knots.at(index)) /
                   (this->knots.at(index + d) - this->knots.at(index)) * *(end - 1);
    }

    return index - this->degree;
  }

  T deboor(size_t index, T value)
  {
    for (size_t j = 0; j <= this->degree; j++)
    {
      this->support[j] = this->control_points.at(j + index - this->degree);
    }

    T alpha = 0;
    for (size_t r = 1; r <= this->degree; r++)
    {
      for (size_t j = this->degree; j >= r; j--)
      {
        alpha = (value - this->knots.at(j + index - this->degree)) /
                (this->knots.at(j + 1 + index - r) - this->knots.at(j + index - this->degree));
        this->support[j] = (1.0 - alpha) * this->support[j - 1] + alpha * this->support[j];
      }
    }

    return this->support[this->degree];
  }
};

} // namespace bsplinex::bspline

#endif
