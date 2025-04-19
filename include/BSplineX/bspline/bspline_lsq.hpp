#ifndef BSPLINEX_BSPLINE_BSPLINE_LSQ_HPP
#define BSPLINEX_BSPLINE_BSPLINE_LSQ_HPP

// Standard includes
#include <functional>
#include <iterator>
#include <vector>

// #ifndef NDEBUG
// #include <iostream>
// #endif

// Third-party includes
#include <Eigen/Dense>
#include <Eigen/Jacobi>
#include <Eigen/src/Core/ConditionEstimator.h>

// For some reason Eigen has a couple of set but unused variables
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-but-set-variable"
#endif

#include <Eigen/Sparse>

#if defined(__clang__)
#pragma clang diagnostic pop
#endif
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

// BSplineX includes
#include "BSplineX/control_points/c_data.hpp"
#include "BSplineX/control_points/control_points.hpp"
#include "BSplineX/defines.hpp"
#include "BSplineX/types.hpp"
#include "BSplineX/views.hpp"
#include "BSplineX/windows.hpp"

namespace bsplinex::lsq
{

using namespace constants;

template <typename T, class MatrixType>
class LSQMatrix
{

  // NOTE: bertolazzi says that the LU algorithm uses roughly half the computations as QR. it is
  // less stable, but for a band matrix it may be fine. plus he suggests to sort the input points
  // as that may improve performance substantially, especially if we develop a specialised LU band
  // algorithm.

public:
  virtual T &operator()(size_t row, size_t col)               = 0;
  virtual Eigen::VectorX<T> solve(Eigen::VectorX<T> const &b) = 0;
  [[nodiscard]] virtual size_t num_rows() const               = 0;
  [[nodiscard]] virtual size_t num_cols() const               = 0;
  [[nodiscard]] virtual T conditioning_number() const         = 0;
};

template <typename T>
class LSQMatrix<T, Eigen::MatrixX<T>>
{
private:
  using mat_t = Eigen::MatrixX<T>;
  using vec_t = Eigen::VectorX<T>;

  mat_t A;

public:
  LSQMatrix(size_t num_rows, size_t num_cols) : A{mat_t::Zero(num_rows, num_cols)} {}

  T &operator()(size_t row, size_t col) { return this->A(row, col); }

  vec_t solve(vec_t const &b)
  {
    Eigen::LeastSquaresConjugateGradient<mat_t> lscg;
    lscg.compute(this->A);
    vec_t x = lscg.solve(b);

    // #ifndef NDEBUG
    //     std::cout << "DENSE" << "\n";
    //     std::cout << "iterations: " << lscg.iterations() << "\n";
    //     std::cout << "error: " << lscg.error() << "\n";
    //     std::cout << "|Ax - b| = " << (this->A * x - b).norm() << "\n";
    //     std::cout << "|b| = " << b.norm() << "\n";
    //     std::cout << "|Ax - b| / |b| = " << (this->A * x - b).norm() / b.norm() << "\n";
    //     std::cout << "condition number = " << this->conditioning_number() << std::endl;
    // #endif
    return x;
  }

  [[nodiscard]] size_t num_rows() const { return this->A.rows(); }

  [[nodiscard]] size_t num_cols() const { return this->A.cols(); }

  [[nodiscard]] T conditioning_number() const
  {
    Eigen::JacobiSVD<mat_t> const svd(this->A);
    T cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);

    return cond;
  }
};

template <typename T>
class LSQMatrix<T, Eigen::SparseMatrix<T>>
{
private:
  using mat_t = Eigen::SparseMatrix<T>;
  using vec_t = Eigen::VectorX<T>;
  mat_t A;

public:
  LSQMatrix(size_t num_rows, size_t num_cols, size_t num_nnz) : A(num_rows, num_cols)
  {
    this->A.reserve(num_nnz);
  }

  T &operator()(size_t row, size_t col) { return this->A.coeffRef(row, col); }

  vec_t solve(vec_t const &b)
  {
    this->A.makeCompressed();

    Eigen::LeastSquaresConjugateGradient<mat_t> lscg;
    lscg.compute(this->A);
    vec_t x = lscg.solve(b);

    // #ifndef NDEBUG
    //     std::cout << "SPARSE\n";
    //     std::cout << "iterations: " << lscg.iterations() << "\n";
    //     std::cout << "error: " << lscg.error() << "\n";
    //     std::cout << "|Ax - b| = " << (this->A * x - b).norm() << "\n";
    //     std::cout << "|b| = " << b.norm() << "\n";
    //     std::cout << "|Ax - b| / |b| = " << (this->A * x - b).norm() / b.norm() << std::endl;
    // #endif

    return x;
  }

  [[nodiscard]] size_t num_rows() const { return A.rows(); }

  [[nodiscard]] size_t num_cols() const { return A.cols(); }

  [[nodiscard]] T conditioning_number() const
  {
    // Eigen::JacobiSVD<Eigen::SparseMatrix<T>> svd(this->A);
    // T cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);

    return 0;
  }
};

template <typename T>
struct Condition
{
  T x_value;
  T y_value;
  size_t derivative_order;

  Condition(T x_value, T y_value, size_t derivative_order)
      : x_value{x_value}, y_value{y_value}, derivative_order{derivative_order}
  {
  }
};

template <typename T, class Iter>
std::vector<Condition<T>> create_sorted_conditions(
    views::ArrayView<Iter> const &x,
    views::ArrayView<Iter> const &y,
    std::vector<Condition<T>> const &additional_conditions
)
{
  debugassert(x.size() == y.size(), "x and y must have the same size.");

  std::vector<Condition<T>> conditions;
  conditions.reserve(x.size() + additional_conditions.size());

  for (size_t k{0}; k < x.size(); k++)
  {
    conditions.emplace_back(x[k], y[k], 0);
  }
  conditions.insert(conditions.end(), additional_conditions.begin(), additional_conditions.end());

  std::sort(
      conditions.begin(),
      conditions.end(),
      [&](Condition<T> const &a, Condition<T> const &b) { return a.x_value < b.x_value; }
  );

  return conditions;
}

template <class LSQMatrix, typename T, class Iter, BoundaryCondition BC>
void fill(
    LSQMatrix &A,
    Eigen::VectorX<T> &b,
    size_t const degree,
    std::function<size_t(T, size_t, std::vector<T> &)> nnz_basis,
    views::ArrayView<Iter> const &x,
    views::ArrayView<Iter> const &y,
    std::vector<Condition<T>> const &additional_conditions
)
{
  debugassert(x.size() == y.size(), "x and y must have the same size.");
  debugassert(
      A.num_rows() == x.size() + additional_conditions.size(),
      "lsq_matrix must have the same size as x + additional_conditions."
  );

  size_t const num_rows = A.num_rows();
  size_t const num_cols = A.num_cols();

  std::vector<Condition<T>> conditions = create_sorted_conditions(x, y, additional_conditions);

  std::vector<T> nnz(degree + 1, ZERO<T>);
  for (size_t i{0}; i < num_rows; i++)
  {
    Condition<T> const &condition = conditions.at(i);

    size_t const index = nnz_basis(condition.x_value, condition.derivative_order, nnz);
    for (size_t j{0}; j <= degree; j++)
    {
      if constexpr (BoundaryCondition::PERIODIC == BC)
      {
        // TODO: avoid modulo
        A(i, (j + index) % num_cols) += nnz.at(j);
      }
      else
      {
        A(i, j + index) = nnz.at(j);
      }
    }
    b(i) = condition.y_value;

    std::fill(nnz.begin(), nnz.end(), ZERO<T>);
  }

  // Check the conditioning number
  // std::cout << "conditioning number for A is: " << A.conditioning_number() << std::endl;
}

template <typename T, class Iter, BoundaryCondition BC>
control_points::ControlPoints<T, BC>
lsq(size_t const degree,
    size_t const knots_size,
    std::function<size_t(T, size_t, std::vector<T> &)> nnz_basis,
    views::ArrayView<Iter> const &x,
    views::ArrayView<Iter> const &y,
    std::vector<Condition<T>> const &additional_conditions)
{
  debugassert(x.size() == y.size(), "x and y must have the same size.");

  using iter_type = typename std::vector<T>::const_iterator;

  size_t const num_rows{x.size() + additional_conditions.size()};
  size_t num_cols{knots_size - degree - 1};

  if constexpr (BoundaryCondition::PERIODIC == BC)
  {
    num_cols -= degree;
  }

  Eigen::VectorX<T> b(num_rows);

  if (num_cols > constants::DENSE_MAX_COL)
  {
    using sparse_lsq = LSQMatrix<T, Eigen::SparseMatrix<T>>;

    sparse_lsq A(num_rows, num_cols, num_cols * (degree + 1));
    fill<sparse_lsq, T, iter_type, BC>(A, b, degree, nnz_basis, x, y, additional_conditions);
    Eigen::VectorX<T> res = A.solve(b);

    using difference_type = typename Eigen::VectorX<T>::iterator::difference_type;

    auto const res_size = static_cast<difference_type>(res.size());
    return control_points::ControlPoints<T, BC>{
        control_points::Data<T>{std::vector<T>{res.data(), std::next(res.data(), res_size)}}, degree
    };
  }
  else
  {
    using dense_lsq = LSQMatrix<T, Eigen::MatrixX<T>>;

    dense_lsq A(num_rows, num_cols);
    fill<dense_lsq, T, iter_type, BC>(A, b, degree, nnz_basis, x, y, additional_conditions);
    Eigen::VectorX<T> res = A.solve(b);

    using difference_type = typename Eigen::VectorX<T>::iterator::difference_type;

    auto const res_size = static_cast<difference_type>(res.size());
    return control_points::ControlPoints<T, BC>{
        control_points::Data<T>{std::vector<T>{res.data(), std::next(res.data(), res_size)}}, degree
    };
  }
}

} // namespace bsplinex::lsq

#endif
