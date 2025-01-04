#ifndef BSPLINE_HPP
#define BSPLINE_HPP

// Standard includes
#include <algorithm>
#include <functional>
#include <sstream>
#include <vector>

// Third-party includes
#include <Eigen/Dense>

// For some reason Eigen has a couple of set but unused variables
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-but-set-variable"
#endif
#include <Eigen/Sparse>
#ifdef __clang__
#pragma clang diagnostic pop
#endif
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

// BSplineX includes
#include "BSplineX/control_points/control_points.hpp"
#include "BSplineX/defines.hpp"
#include "BSplineX/knots/knots.hpp"
#include "BSplineX/types.hpp"

namespace bsplinex::bspline
{

template <typename T, Curve C, BoundaryCondition BC, Extrapolation EXT>
class BSpline
{
private:
  knots::Knots<T, C, BC, EXT> knots{};
  control_points::ControlPoints<T, BC> control_points{};
  size_t degree{};
  std::vector<T> support{};

public:
  BSpline() { DEBUG_LOG_CALL(); }

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

  BSpline(BSpline const &other)
      : knots(other.knots), control_points(other.control_points), degree(other.degree),
        support(other.support)
  {
    DEBUG_LOG_CALL();
  }

  BSpline(BSpline &&other) noexcept
      : knots(std::move(other.knots)), control_points(std::move(other.control_points)),
        degree(other.degree), support(std::move(other.support))
  {
    DEBUG_LOG_CALL();
  }

  ~BSpline() noexcept { DEBUG_LOG_CALL(); }

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

  T evaluate(T value)
  {
    auto index_value_pair = this->knots.find(value);
    return this->deboor(index_value_pair.first, index_value_pair.second);
  }

  std::vector<T> basis(T value)
  {
    std::vector<T> basis_functions(this->degree + 1, (T)0);

    size_t index = this->nnz_basis(
        this->degree, this->knots, value, basis_functions.begin(), basis_functions.end()
    );

    basis_functions.insert(basis_functions.begin(), index, (T)0);
    basis_functions.insert(
        basis_functions.end(), this->control_points.size() - index - this->degree - 1, (T)0
    );

    return basis_functions;
  }

  template <typename It>
  static size_t nnz_basis(
      size_t degree,
      knots::Knots<T, C, BC, EXT> const &knots,
      T value,
      [[maybe_unused]] It begin,
      It end
  )
  {
    assertm(
        (end - begin) == (long long)(degree + 1),
        "Unexpected number of basis asked, exactly degree + 1 basis can be asked"
    );

    assertm(
        std::all_of(begin, end, [](T i) { return (T)0 == i; }),
        "Initial basis must be initialised to zero"
    );

    auto [index, val] = knots.find(value);

    *(end - 1) = 1.0;
    for (size_t d{1}; d <= degree; d++)
    {
      *(end - 1 - d) = (knots.at(index + 1) - val) /
                       (knots.at(index + 1) - knots.at(index - d + 1)) * *(end - 1 - d + 1);
      for (size_t i{index - d + 1}; i < index; i++)
      {
        *(end - 1 - index + i) =
            (val - knots.at(i)) / (knots.at(i + d) - knots.at(i)) * *(end - 1 - index + i) +
            (knots.at(i + d + 1) - val) / (knots.at(i + d + 1) - knots.at(i + 1)) *
                *(end - 1 - index + i + 1);
      }
      *(end - 1) = (val - knots.at(index)) / (knots.at(index + d) - knots.at(index)) * *(end - 1);
    }

    return index - degree;
  }

  void fit(std::vector<T> const &x, std::vector<T> const &y)
  {
    if (x.size() != y.size())
    {
      throw std::runtime_error("x and y must have the same size");
    }

    // NOTE: bertolazzi says that the LU algorithm uses roughly half the computations as QR. it is
    // less stable, but for a band matrix it may be fine. plus he suggests to sort the input points
    // as that may improve performance substantially, especially if we develop a specialised LU band
    // algorithm.

    std::vector<T> nnz_basis_vec(this->degree + 1, (T)0);
    Eigen::Map<Eigen::VectorX<T> const> b(y.data(), y.size());
    Eigen::VectorX<T> res;
    size_t num_cols{this->control_points.size()};
    if constexpr (BoundaryCondition::PERIODIC == BC)
    {
      num_cols -= this->degree;
    }

    if (num_cols <= DENSE_MAX_COL)
    {
      Eigen::MatrixX<T> A = Eigen::MatrixX<T>::Zero(x.size(), num_cols);

      size_t index{0};
      for (size_t i{0}; i < x.size(); i++)
      {
        index = this->nnz_basis(
            this->degree, this->knots, x.at(i), nnz_basis_vec.begin(), nnz_basis_vec.end()
        );
        for (size_t j{0}; j <= this->degree; j++)
        {
          // TODO: avoid modulo
          A(i, (j + index) % num_cols) += nnz_basis_vec.at(j);
        }
        std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
      }

      res = A.colPivHouseholderQr().solve(b);
    }
    else
    {
      Eigen::SparseMatrix<T> A(x.size(), num_cols);
      A.reserve(num_cols * (this->degree + 1));

      size_t index{0};
      for (size_t i{0}; i < x.size(); i++)
      {
        index = this->nnz_basis(
            this->degree, this->knots, x.at(i), nnz_basis_vec.begin(), nnz_basis_vec.end()
        );
        for (size_t j{0}; j <= this->degree; j++)
        {
          A.coeffRef(i, (j + index) % num_cols) += nnz_basis_vec.at(j);
        }
        std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
      }
      A.makeCompressed();

      Eigen::SparseQR<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>> solver{};
      solver.compute(A);
      res = solver.solve(b);
    }

    new (&this->control_points) control_points::ControlPoints<T, BC>{
        {{res.data(), res.data() + res.rows() * res.cols()}}, this->degree
    };
  }

  void interpolate(
      std::vector<T> const &x,
      std::vector<T> const &y,
      [[maybe_unused]] std::vector<T> const &padding
  )
  {
    if constexpr (Curve::UNIFORM == C)
    {
      return;
    }
    else
    {
      if (x.size() != y.size())
      {
        throw std::runtime_error("x and y must have the same size");
      }
      // NOTE: thank the STL for this wonderful backwards built sort check. Think it as if std::less
      // is <= and std::less_equal is <.
      if (!std::is_sorted(x.begin(), x.end(), std::less_equal<T>{}))
      {
        throw std::runtime_error("x must be sorted w.r.t. operator <");
      }

      if constexpr (BoundaryCondition::OPEN == BC)
      {
        if (padding.size() != 2 * this->degree)
        {
          throw std::runtime_error("padding must be = 2 * degree");
        }
        size_t open_knots_size = x.size() + 2 * this->degree;
        std::vector<T> open_knots{};
        open_knots.reserve(open_knots_size);
        for (size_t i{0}; i < this->degree; i++)
        {
          open_knots.push_back(padding.at(i));
        }
        for (auto const &elem : x)
        {
          open_knots.push_back(elem);
        }
        for (size_t i{0}; i < this->degree; i++)
        {
          open_knots.push_back(padding.at(i + this->degree));
        }

        new (&this->knots) knots::Knots<T, C, BC, EXT>{{open_knots}, this->degree};
      }
      else
      {
        new (&this->knots) knots::Knots<T, C, BC, EXT>{{x}, this->degree};
      }

      size_t num_cols = x.size() + this->degree - 1;
      size_t eq_aft   = (this->degree - 1) / 2;
      size_t eq_bef   = (this->degree - 1) - eq_aft;
      if constexpr (BoundaryCondition::PERIODIC == BC)
      {
        num_cols -= this->degree - 1;
        eq_aft    = 0;
        eq_bef    = 0;
      }
      std::vector<T> nnz_basis_vec(this->degree + 1, (T)0);
      Eigen::MatrixX<T> A = Eigen::MatrixX<T>::Zero(num_cols, num_cols);
      Eigen::VectorX<T> b(num_cols);
      Eigen::VectorX<T> res;

      size_t index{0};
      size_t i{0};

      for (; i < eq_bef; i++)
      {
        A(i, i) += (T)1;
        b(i)     = y.front();
      }

      for (; i < num_cols - eq_aft; i++)
      {
        index = this->nnz_basis(
            this->degree, this->knots, x.at(i - eq_bef), nnz_basis_vec.begin(), nnz_basis_vec.end()
        );
        for (size_t j{0}; j <= this->degree; j++)
        {
          A(i, (j + index) % num_cols) += nnz_basis_vec.at(j);
        }
        std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);

        b(i) = y.at(i - eq_bef);
      }

      for (; i < num_cols; i++)
      {
        A(i, i) += (T)1;
        b(i)     = y.back();
      }

      res = A.colPivHouseholderQr().solve(b);

      new (&this->control_points) control_points::ControlPoints<T, BC>{
          {{res.data(), res.data() + res.rows() * res.cols()}}, this->degree
      };
    }
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

    throw std::runtime_error(ss.str());
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
