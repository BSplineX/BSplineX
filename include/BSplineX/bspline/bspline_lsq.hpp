#ifndef BSPLINE_LSQ_HPP
#define BSPLINE_LSQ_HPP

// Standard includes
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
#include "BSplineX/bspline/bspline.hpp"
#include "BSplineX/control_points/control_points.hpp"
#include "BSplineX/defines.hpp"
#include "BSplineX/knots/knots.hpp"
#include "BSplineX/types.hpp"

namespace bsplinex::lsq
{

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
  virtual size_t num_rows() const                             = 0;
  virtual size_t num_cols() const                             = 0;
};

template <typename T>
class LSQMatrix<T, Eigen::MatrixX<T>>
{
private:
  Eigen::MatrixX<T> A;

public:
  LSQMatrix(size_t num_rows, size_t num_cols) : A{Eigen::MatrixX<T>::Zero(num_rows, num_cols)} {}

  T &operator()(size_t row, size_t col) { return A(row, col); }

  Eigen::VectorX<T> solve(Eigen::VectorX<T> const &b) { return A.colPivHouseholderQr().solve(b); }

  size_t num_rows() const { return A.num_rows(); }

  size_t num_cols() const { return A.num_cols(); }
};

template <typename T>
class LSQMatrix<T, Eigen::SparseMatrix<T>>
{
private:
  Eigen::SparseMatrix<T> A;

public:
  LSQMatrix(size_t num_rows, size_t num_cols, size_t num_nnz) : A(num_rows, num_cols)
  {
    A.reserve(num_nnz);
  }

  T &operator()(size_t row, size_t col) { return A.coeffRef(row, col); }

  Eigen::VectorX<T> solve(Eigen::VectorX<T> const &b)
  {
    A.makeCompressed();

    Eigen::SparseQR<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>> solver{};
    solver.compute(A);
    return solver.solve(b);
  }

  size_t num_rows() const { return A.num_rows(); }

  size_t num_cols() const { return A.num_cols(); }
};

template <typename T>
struct Condition
{
  T x_value;
  T y_value;
  size_t derivative_degree;

  Condition(T x_value, T y_value, size_t derivative_degree)
      : x_value{x_value}, y_value{y_value}, derivative_degree{derivative_degree}
  {
  }
};

template <class LSQMatrix, typename T, Curve C, BoundaryCondition BC, Extrapolation EXT>
void fill(
    LSQMatrix &A,
    Eigen::VectorX<T> &b,
    size_t degree,
    knots::Knots<T, C, BC, EXT> const &knots,
    std::vector<T> const &x,
    std::vector<T> const &y,
    std::vector<Condition<T>> const &additional_conditions
)
{
  assertm(x.size() == y.size(), "x and y must have the same size.");
  assertm(
      A.num_rows() == x.size() + additional_conditions.size(),
      "lsq_matrix must have the same size as x + additional_conditions."
  );

  size_t const &num_rows = A.num_rows();
  size_t const &num_cols = A.num_cols();

  std::vector<Condition<T>> conditions;
  conditions.reserve(num_rows);

  size_t k{0};
  for (; k < x.size(); k++)
  {
    conditions.emplace_back(x[k], y[k], 0);
  }
  for (k -= x.size(); k < additional_conditions.size(); k++)
  {
    conditions.push_back(additional_conditions[k]);
  }

  std::sort(
      conditions.begin(),
      conditions.end(),
      [&](Condition<T> const &a, Condition<T> const &b) { return a.x_value < b.x_value; }
  );

  std::vector<T> nnz(degree + 1);
  for (size_t i{0}; i < num_rows; i++)
  {
    Condition<T> const &condition = conditions.at(i);

    size_t index = bspline::BSpline<T, C, BC, EXT>::nnz_basis(
        degree, knots, condition.x_value, nnz.begin(), nnz.end()
    );
    for (size_t j{0}; j <= degree; j++)
    {
      if constexpr (BoundaryCondition::PERIODIC == BC)
      {
        // TODO: avoid modulo
        A.coeffRef(i, (j + index) % num_cols) += nnz.at(j);
      }
      else
      {
        A.coeffRef(i, j + index) = nnz.at(j);
      }
    }
    b(i) = condition.y_value;

    std::fill(nnz.begin(), nnz.end(), (T)0);
  }
}

template <typename T, Curve C, BoundaryCondition BC, Extrapolation EXT>
control_points::ControlPoints<T, BC>
lsq(size_t degree,
    knots::Knots<T, C, BC, EXT> const &knots,
    std::vector<T> const &x,
    std::vector<T> const &y,
    std::vector<Condition<T>> const &additional_conditions)
{
  if (x.size() != y.size())
  {
    throw std::runtime_error("x and y must have the same size");
  }

  size_t const &num_rows = x.size();
  size_t num_cols{knots.size() - degree - 1};
  if constexpr (BoundaryCondition::PERIODIC == BC)
  {
    num_cols -= degree;
  }

  Eigen::VectorX<T> b;
  b.reserve(num_rows);

  if (num_cols > DENSE_MAX_COL)
  {
    LSQMatrix<T, Eigen::SparseMatrix<T>> A(num_rows, num_cols, num_cols * (degree + 1));
    fill(A, b, degree, knots, x, y, additional_conditions);
    Eigen::VectorX<T> res = A.solve(b);

    return control_points::ControlPoints<T, BC>{
        {{res.data(), res.data() + res.rows() * res.cols()}}, degree
    };
  }
  else
  {
    LSQMatrix<T, Eigen::MatrixX<T>> A(num_rows, num_cols);
    fill(A, b, degree, knots, x, y, additional_conditions);
    Eigen::VectorX<T> res = A.solve(b);

    return control_points::ControlPoints<T, BC>{
        {{res.data(), res.data() + res.rows() * res.cols()}}, degree
    };
  }
}

} // namespace bsplinex::lsq

#endif
