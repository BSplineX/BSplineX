#ifndef BSPLINEX_BSPLINE_BSPLINE_HPP
#define BSPLINEX_BSPLINE_BSPLINE_HPP

// Standard includes
#include <algorithm>
#include <limits>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>

// BSplineX includes
#include "BSplineX/bspline/bspline_lsq.hpp"
#include "BSplineX/control_points/control_points.hpp"
#include "BSplineX/defines.hpp"
#include "BSplineX/knots/knots.hpp"
#include "BSplineX/types.hpp"
#include "BSplineX/views.hpp"

namespace bsplinex::bspline
{

using namespace constants;

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
public:
  static Curve const curve_type{C};

  static BoundaryCondition const boundary_condition_type{BC};

  static Extrapolation const extrapolation_type{EXT};

private:
  using vec_const_iter = typename std::vector<T>::const_iterator;
  using vec_const_view = typename views::ArrayView<vec_const_iter>;
  using vec_iter       = typename std::vector<T>::iterator;
  using vec_view       = typename views::ArrayView<vec_iter>;

  knots::Knots<T, C, BC, EXT> knots{};
  control_points::ControlPoints<T, BC> control_points{};
  size_t degree{};
  std::vector<T> support{};
  std::unique_ptr<BSpline> derivative;

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
   * @param derivative_order The order of the derivative to evaluate (0 is the function itself).
   * @return The value of the BSpline at the given value.
   */
  T evaluate(T value, size_t derivative_order = 0)
  {
    releaseassert(derivative_order <= this->degree, "derivative_order must be <= degree");
    if (0 == derivative_order)
    {
      auto [index, x_value] = this->knots.find(value);
      return this->deboor(index, x_value);
    }
    else
    {
      if (not this->derivative)
      {
        this->compute_derivative();
      }

      return derivative->evaluate(value, derivative_order - 1);
    }
  }

  /**
   * @brief Evaluates the BSpline at the given values.
   *
   * @param values The values to evaluate the BSpline at.
   * @param derivative_order The order of the derivative to evaluate (0 is the function itself).
   * @return The values of the BSpline at the given values.
   */
  std::vector<T> evaluate(std::vector<T> const &values, size_t derivative_order = 0)
  {
    std::vector<T> results;
    results.reserve(values.size());

    for (auto const &value : values)
    {
      results.push_back(this->evaluate(value, derivative_order));
    }

    return results;
  }

  /**
   * @brief Computes the basis functions at a given value.
   *
   * @param value The value to compute the basis functions at.
   * @return The basis functions at the given value.
   */
  std::vector<T> basis(T value)
  {
    auto [index, basis_functions] = this->nnz_basis(value, 0);

    basis_functions.insert(basis_functions.begin(), index, ZERO<T>);
    basis_functions.insert(
        basis_functions.end(), this->control_points.size() - index - this->degree - 1, ZERO<T>
    );

    return basis_functions;
  }

  std::pair<size_t, std::vector<T>> nnz_basis(T value, size_t derivative_order)
  {
    std::vector<T> nnz(this->degree + 1, ZERO<T>);
    size_t index = this->nnz_basis<vec_iter>(value, derivative_order, {nnz.begin(), nnz.end()});

    return std::make_pair(index, std::move(nnz));
  }

  /**
   * @brief Gives the (inclusive) boundary of the BSpline domain (i.e., [x_min, x_max])
   */
  std::pair<T, T> domain() const { return this->knots.domain(); }

  /**
   * @brief Fits the BSpline to some data.
   *
   * @param x The x values of the data.
   * @param y The y values of the data.
   */
  void fit(std::vector<T> const &x, std::vector<T> const &y)
  {
    releaseassert(x.size() == y.size(), "x and y must have the same size");

    this->control_points = std::move(
        lsq::lsq<T, vec_const_iter, BC>(
            degree,
            knots.size(),
            [this](T value, size_t derivative_order, std::vector<T> &vec) -> size_t
            {
              return this->nnz_basis<vec_iter>(value, derivative_order, {vec.begin(), vec.end()});
            },
            vec_const_view{x.begin(), x.end()},
            vec_const_view{y.begin(), y.end()},
            {}
        )
    );
    this->invalidate_derivative();
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
    this->invalidate_derivative();

    knots::Knots<T, C, BC, EXT> new_knots;
    if constexpr (Curve::UNIFORM == C)
    {
      T step = x.at(1) - x.at(0);
      for (size_t i{0}; i < x.size() - 1; i++)
      {
        releaseassert(
            std::abs(x.at(i + 1) - x.at(i) - step) <=
                std::numeric_limits<T>::epsilon() * 1000.0, // HACK:
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

    vec_const_view x_view, y_view;
    if constexpr (BoundaryCondition::OPEN == BC)
    {
      x_view = vec_const_view{std::next(x.begin(), degree), std::prev(x.end(), degree)};
      y_view = vec_const_view{std::next(y.begin(), degree), std::prev(y.end(), degree)};
    }
    else if constexpr (BoundaryCondition::CLAMPED == BC)
    {
      x_view = vec_const_view{x.begin(), x.end()};
      y_view = vec_const_view{y.begin(), y.end()};
    }
    else if constexpr (BoundaryCondition::PERIODIC == BC)
    {
      x_view = vec_const_view{x.begin(), std::prev(x.end(), 1)};
      y_view = vec_const_view{y.begin(), std::prev(y.end(), 1)};
    }
    else
    {
      static_assert(false, "Unknown boundary condition, you should never get here!");
    }

    this->control_points = std::move(
        lsq::lsq<T, vec_const_iter, BC>(
            this->degree,
            this->knots.size(),
            [this](T value, size_t derivative_order, std::vector<T> &vec) -> size_t
            {
              return this->nnz_basis<vec_iter>(value, derivative_order, {vec.begin(), vec.end()});
            },
            x_view,
            y_view,
            additional_conditions
        )
    );
    this->invalidate_derivative();
  }

  std::vector<T> get_control_points() const
  {
    std::vector<T> ctrl_pts{};
    ctrl_pts.reserve(this->control_points.size());

    for (size_t i{0}; i < this->control_points.size(); i++)
    {
      ctrl_pts.push_back(this->control_points.at(i));
    }

    return ctrl_pts;
  }

  std::vector<T> get_knots() const
  {
    std::vector<T> knots{};
    knots.reserve(this->knots.size());
    for (size_t i{0}; i < this->knots.size(); i++)
    {
      knots.push_back(this->knots.at(i));
    }

    return knots;
  }

private:
  BSpline(
      knots::Knots<T, C, BC, EXT> const &knots,
      control_points::ControlPoints<T, BC> const &control_points,
      size_t degree
  )
      : knots{knots}, control_points{control_points}, degree{degree}
  {
    DEBUG_LOG_CALL();
    this->check_sizes();
    this->support.resize(this->degree + 1);
  }

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

  // Algorithm adapted from:
  // - https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve-coef.html
  // - https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-derv.html
  template <typename Iter>
  size_t nnz_basis(T value, size_t derivative_order, vec_view nnz)
  {
    debugassert(
        nnz.size() == this->degree + 1,
        "Unexpected number of basis asked, exactly degree + 1 basis can be asked"
    );

    debugassert(
        std::all_of(nnz.begin(), nnz.end(), [](T i) { return ZERO<T> == i; }),
        "Initial basis must be initialised to zero"
    );

    releaseassert(this->degree > derivative_order, "Asked for derivative_order >= degree");

    auto [index, val] = this->knots.find(value);

    // Compute the basis of degree - derivative_order
    nnz.back() = 1.0;
    for (size_t d{1}; d <= this->degree - derivative_order; d++)
    {
      size_t const idx{this->degree - d};
      nnz.at(idx) = (this->knots.at(index + 1) - val) /
                    (this->knots.at(index + 1) - this->knots.at(index - d + 1)) * nnz.at(idx + 1);
      for (size_t i{index - d + 1}; i < index; ++i)
      {
        size_t const idx_in{this->degree - index};

        T const den_1{(this->knots.at(i + d) - this->knots.at(i))};
        T const den_2{(this->knots.at(i + d + 1) - this->knots.at(i + 1))};
        T const num_1{val - this->knots.at(i)};
        T const num_2{this->knots.at(i + d + 1) - val};
        T const basis_1{num_1 / den_1 * nnz.at(idx_in + i)};
        T const basis_2{num_2 / den_2 * nnz.at(idx_in + i + 1)};

        nnz.at(idx_in + i) = basis_1 + basis_2;
      }

      T const den_1{this->knots.at(index + d) - this->knots.at(index)};
      T const num_1{val - this->knots.at(index)};

      nnz.back() = num_1 / den_1 * nnz.back();
    }

    // Compute the derivatives up to derivative_order
    for (size_t p{this->degree + 1 - derivative_order}; p <= this->degree; p++)
    {
      for (size_t i{0}; i < this->degree; i++)
      {
        size_t const idx{index - this->degree + i};

        T const den_1{this->knots.at(idx + p) - this->knots.at(idx)};
        T const den_2{this->knots.at(idx + p + 1) - this->knots.at(idx + 1)};
        T const base_1 = den_1 ? p / den_1 * (nnz.at(i)) : ZERO<T>;
        T const base_2 = den_2 ? p / den_2 * (nnz.at(i + 1)) : ZERO<T>;

        nnz.at(i) = base_1 - base_2;
      }

      T const den_1{this->knots.at(index + p) - this->knots.at(index)};

      nnz.back() = den_1 ? p / den_1 * nnz.back() : ZERO<T>;
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
      for (size_t j = this->degree; j >= r; --j)
      {
        alpha = (value - this->knots.at(j + index - this->degree)) /
                (this->knots.at(j + 1 + index - r) - this->knots.at(j + index - this->degree));
        this->support[j] = (1.0 - alpha) * this->support[j - 1] + alpha * this->support[j];
      }
    }

    return this->support[this->degree];
  }

  void compute_derivative()
  {
    debugassert(this->degree > 0, "Cannot compute derivative of a 0-degree bspline");
    this->derivative = std::unique_ptr<BSpline>(new BSpline(
        this->knots.get_derivative_knots(),
        this->control_points.get_derivative_control_points(this->knots),
        this->degree - 1
    ));
  }

  void invalidate_derivative() { this->derivative = nullptr; }
};

} // namespace bsplinex::bspline

#endif
