#ifndef BSPLINE_FIT_HPP
#define BSPLINE_FIT_HPP

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

namespace bsplinex::fit
{

template <typename T, Curve C, BoundaryCondition BC, Extrapolation EXT>
control_points::ControlPoints<T, BC>
fit(size_t degree,
    knots::Knots<T, C, BC, EXT> const &knots,
    std::vector<T> const &x,
    std::vector<T> const &y)
{
  if ((knots.size() - degree - 1) > DENSE_MAX_COL)
  {
    return sparse(degree, knots, x, y);
  }

  return dense(degree, knots, x, y);
}

// NOTE: bertolazzi says that the LU algorithm uses roughly half the computations as QR. it is
// less stable, but for a band matrix it may be fine. plus he suggests to sort the input points
// as that may improve performance substantially, especially if we develop a specialised LU band
// algorithm.

template <typename T, Curve C, BoundaryCondition BC, Extrapolation EXT>
control_points::ControlPoints<T, BC> dense(
    size_t degree,
    knots::Knots<T, C, BC, EXT> const &knots,
    std::vector<T> const &x,
    std::vector<T> const &y
)
{
  if (x.size() != y.size())
  {
    throw std::runtime_error("x and y must have the same size");
  }

  size_t num_cols{knots.size() - degree - 1};
  if constexpr (BoundaryCondition::PERIODIC == BC)
  {
    num_cols -= degree;
  }

  size_t num_rows     = x.size();
  Eigen::MatrixX<T> A = Eigen::MatrixX<T>::Zero(num_rows, num_cols);
  Eigen::VectorX<T> b;
  b.reserve(num_rows);

  size_t index{0};
  std::vector<T> nnz_basis_vec(degree + 1, (T)0);
  T initial_x = x.front();
  T final_x   = x.back();

  for (size_t i{0}; i < num_rows; i++)
  {
    index = bspline::BSpline<T, C, BC, EXT>::nnz_basis(
        degree, knots, x.at(i), nnz_basis_vec.begin(), nnz_basis_vec.end()
    );
    for (size_t j{0}; j <= degree; j++)
    {
      if constexpr (BoundaryCondition::PERIODIC == BC)
      {
        // TODO: avoid modulo
        A(i, (j + index) % num_cols) += nnz_basis_vec.at(j);
      }
      else
      {
        A(i, j + index) = nnz_basis_vec.at(j);
      }
    }
    b(i) = y.at(i);

    std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
  }

  Eigen::VectorX<T> res = A.colPivHouseholderQr().solve(b);

  return control_points::ControlPoints<T, BC>{
      {{res.data(), res.data() + res.rows() * res.cols()}}, degree
  };
}

template <typename T, Curve C, BoundaryCondition BC, Extrapolation EXT>
control_points::ControlPoints<T, BC> sparse(
    size_t degree,
    knots::Knots<T, C, BC, EXT> const &knots,
    std::vector<T> const &x,
    std::vector<T> const &y
)
{
  if (x.size() != y.size())
  {
    throw std::runtime_error("x and y must have the same size");
  }

  size_t num_cols{knots.size() - degree - 1};
  if constexpr (BoundaryCondition::PERIODIC == BC)
  {
    num_cols -= degree;
  }

  size_t num_rows = x.size();
  Eigen::SparseMatrix<T> A(num_rows, num_cols);
  A.reserve(num_cols * (degree + 1));
  Eigen::VectorX<T> b;
  b.reserve(num_rows);

  size_t index{0};
  std::vector<T> nnz_basis_vec(degree + 1, (T)0);
  T initial_x = x.front();
  T final_x   = x.back();

  for (size_t i{0}; i < x.size(); i++)
  {
    index = bspline::BSpline<T, C, BC, EXT>::nnz_basis(
        degree, knots, x.at(i), nnz_basis_vec.begin(), nnz_basis_vec.end()
    );
    for (size_t j{0}; j <= degree; j++)
    {
      if constexpr (BoundaryCondition::PERIODIC == BC)
      {
        // TODO: avoid modulo
        A.coeffRef(i, (j + index) % num_cols) += nnz_basis_vec.at(j);
      }
      else
      {
        A.coeffRef(i, j + index) = nnz_basis_vec.at(j);
      }
    }
    b(i) = y.at(i);

    std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
  }
  A.makeCompressed();

  Eigen::SparseQR<Eigen::SparseMatrix<T>, Eigen::COLAMDOrdering<int>> solver{};
  solver.compute(A);
  Eigen::VectorX<T> res = solver.solve(b);

  return control_points::ControlPoints<T, BC>{
      {{res.data(), res.data() + res.rows() * res.cols()}}, degree
  };
}

} // namespace bsplinex::fit

#endif
