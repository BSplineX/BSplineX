#ifndef BSPLINE_INTERPOLATE_HPP
#define BSPLINE_INTERPOLATE_HPP

// Standard includes
#include <algorithm>
#include <functional>
#include <utility>
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

namespace bsplinex::interpolate
{

template <typename T, BoundaryCondition BC, Extrapolation EXT>
std::pair<knots::Knots<T, Curve::UNIFORM, BC, EXT>, control_points::ControlPoints<T, BC>> uniform(
    size_t degree,
    T begin,
    T end,
    size_t num_elems,
    std::vector<T> const &y,
    std::vector<std::pair<size_t, T>> const &initial_conditions,
    std::vector<std::pair<size_t, T>> const &final_conditions
)
{
  size_t step_size = (end - begin) / (num_elems - 1);
  knots::Knots<T, Curve::UNIFORM, BC, EXT> knots{};
  if constexpr (BoundaryCondition::OPEN == BC)
  {
    new (&knots) knots::Knots<T, Curve::UNIFORM, BC, EXT>{
        {begin - degree * step_size, end + degree * step_size, num_elems + 2 * degree}, degree
    };
  }
  else
  {
    new (&knots) knots::Knots<T, Curve::UNIFORM, BC, EXT>{{begin, end, num_elems}, degree};
  }

  control_points::ControlPoints<T, BC> ctrl_pts{};

  if ((knots.size() - degree - 1) > DENSE_MAX_COL)
  {
    ctrl_pts = sparse(degree, knots, y, initial_conditions, final_conditions);
  }
  else
  {
    ctrl_pts = dense(degree, knots, y, initial_conditions, final_conditions);
  }

  return std::make_pair(std::move(knots), std::move(ctrl_pts));
}

template <typename T, BoundaryCondition BC, Extrapolation EXT>
std::pair<knots::Knots<T, Curve::NON_UNIFORM, BC, EXT>, control_points::ControlPoints<T, BC>>
non_uniform(
    size_t degree,
    std::vector<T> const &x,
    std::vector<T> const &y,
    [[maybe_unused]] std::vector<T> const &padding,
    std::vector<std::pair<size_t, T>> const &initial_conditions,
    std::vector<std::pair<size_t, T>> const &final_conditions
)
{
  // NOTE: thank the STL for this wonderful backwards built sort check. Think it as if std::less
  // is <= and std::less_equal is <.
  if (!std::is_sorted(x.begin(), x.end(), std::less_equal<T>{}))
  {
    throw std::runtime_error("x must be sorted w.r.t. operator <");
  }

  knots::Knots<T, Curve::NON_UNIFORM, BC, EXT> knots{};
  if constexpr (BoundaryCondition::OPEN == BC)
  {
    if (padding.size() != 2 * degree)
    {
      throw std::runtime_error("padding must be = 2 * degree");
    }
    size_t open_knots_size = x.size() + 2 * degree;
    std::vector<T> open_knots{};
    open_knots.reserve(open_knots_size);
    for (size_t i{0}; i < degree; i++)
    {
      open_knots.push_back(padding.at(i));
    }
    for (auto const &elem : x)
    {
      open_knots.push_back(elem);
    }
    for (size_t i{0}; i < degree; i++)
    {
      open_knots.push_back(padding.at(i + degree));
    }

    new (&knots) knots::Knots<T, Curve::NON_UNIFORM, BC, EXT>{{open_knots}, degree};
  }
  else
  {
    new (&knots) knots::Knots<T, Curve::NON_UNIFORM, BC, EXT>{{x}, degree};
  }

  control_points::ControlPoints<T, BC> ctrl_pts{};

  if ((knots.size() - degree - 1) > DENSE_MAX_COL)
  {
    ctrl_pts = sparse(degree, knots, x, y, initial_conditions, final_conditions);
  }
  else
  {
    ctrl_pts = dense(degree, knots, x, y, initial_conditions, final_conditions);
  }

  return std::make_pair(std::move(knots), std::move(ctrl_pts));
}

// NOTE: bertolazzi says that the LU algorithm uses roughly half the computations as QR. it is
// less stable, but for a band matrix it may be fine. plus he suggests to sort the input points
// as that may improve performance substantially, especially if we develop a specialised LU band
// algorithm.

template <typename T, BoundaryCondition BC, Extrapolation EXT>
control_points::ControlPoints<T, BC> dense_uniform(
    size_t degree,
    knots::Knots<T, Curve::UNIFORM, BC, EXT> const &knots,
    std::vector<T> const &y,
    std::vector<std::pair<size_t, T>> const &initial_conditions,
    std::vector<std::pair<size_t, T>> const &final_conditions
)
{
  size_t num_cols{knots.size() - degree - 1};
  if constexpr (BoundaryCondition::PERIODIC == BC)
  {
    num_cols -= degree;
  }

  size_t num_initial = initial_conditions.size();
  size_t num_final   = final_conditions.size();
  size_t num_rows    = y.size() + num_initial + num_final;

  if (num_rows != num_cols)
  {
    throw std::runtime_error(
        "num_rows has to be equal to num_cols (i.e., the system has to be properly determined)"
    );
  }

  Eigen::MatrixX<T> A = Eigen::MatrixX<T>::Zero(num_rows, num_cols);
  Eigen::VectorX<T> b;
  b.reserve(num_rows);

  size_t index{0}, i{0};
  std::vector<T> nnz_basis_vec(degree + 1, (T)0);
  T initial_x = knots.at(degree);
  T final_x   = knots.at(knots.size() - degree);
  for (; i < num_initial; i++)
  {
    // TODO: add degree of derivative, same goes for all ::nnz_basis in this function.
    // TODO: no need to call this method, for interpolation we have no reason to call knots.find
    // since we know the point lands on a knot and that all of the point correspond to a knot. Same
    // goes for all the ::nnz_basis in this function. Moreover, it is actually needed to call it
    // once, since the values will stay the same, only the index changes (but simply follows the
    // points we are iterating over since A is a band matrix, aside from the initial/final
    // conditions).
    index = bspline::BSpline<T, Curve::UNIFORM, BC, EXT>::nnz_basis(
        degree, knots, initial_x, nnz_basis_vec.begin(), nnz_basis_vec.end()
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
    b(i) = initial_conditions.at(i).second;

    std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
  }

  for (; i < y.size(); i++)
  {
    index = bspline::BSpline<T, Curve::UNIFORM, BC, EXT>::nnz_basis(
        degree,
        knots,
        knots.at(i - num_initial + degree),
        nnz_basis_vec.begin(),
        nnz_basis_vec.end()
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
    b(i) = y.at(i - num_initial);

    std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
  }

  for (; i < num_rows; i++)
  {
    index = bspline::BSpline<T, Curve::UNIFORM, BC, EXT>::nnz_basis(
        degree, knots, final_x, nnz_basis_vec.begin(), nnz_basis_vec.end()
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
    b(i) = final_conditions.at(i - (num_rows - num_final)).second;

    std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
  }

  Eigen::VectorX<T> res = A.colPivHouseholderQr().solve(b);

  return control_points::ControlPoints<T, BC>{
      {{res.data(), res.data() + res.rows() * res.cols()}}, degree
  };
}

template <typename T, BoundaryCondition BC, Extrapolation EXT>
control_points::ControlPoints<T, BC> sparse_uniform(
    size_t degree,
    knots::Knots<T, Curve::UNIFORM, BC, EXT> const &knots,
    std::vector<T> const &y,
    std::initializer_list<std::pair<size_t, T>> const &initial_conditions,
    std::initializer_list<std::pair<size_t, T>> const &final_conditions
)
{
  size_t num_cols{knots.size() - degree - 1};
  if constexpr (BoundaryCondition::PERIODIC == BC)
  {
    num_cols -= degree;
  }

  size_t num_initial = initial_conditions.size();
  size_t num_final   = final_conditions.size();
  size_t num_rows    = y.size() + num_initial + num_final;

  if (num_rows != num_cols)
  {
    throw std::runtime_error(
        "num_rows has to be equal to num_cols (i.e., the system has to be properly determined)"
    );
  }

  Eigen::SparseMatrix<T> A(num_rows, num_cols);
  A.reserve(num_cols * (degree + 1));
  Eigen::VectorX<T> b;
  b.reserve(num_rows);

  size_t index{0}, i{0};
  std::vector<T> nnz_basis_vec(degree + 1, (T)0);
  T initial_x = knots.at(degree);
  T final_x   = knots.at(knots.size() - degree);
  for (; i < num_initial; i++)
  {
    // TODO: add degree of derivative, same goes for all ::nnz_basis in this function.
    // TODO: same comment about the no need for knots.find.
    index = bspline::BSpline<T, Curve::UNIFORM, BC, EXT>::nnz_basis(
        degree, knots, initial_x, nnz_basis_vec.begin(), nnz_basis_vec.end()
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
    b(i) = initial_conditions.at(i).second;

    std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
  }

  for (; i < y.size(); i++)
  {
    index = bspline::BSpline<T, Curve::UNIFORM, BC, EXT>::nnz_basis(
        degree,
        knots,
        knots.at(i - num_initial + degree),
        nnz_basis_vec.begin(),
        nnz_basis_vec.end()
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
    b(i) = y.at(i - num_initial);

    std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
  }

  for (; i < num_rows; i++)
  {
    index = bspline::BSpline<T, Curve::UNIFORM, BC, EXT>::nnz_basis(
        degree, knots, final_x, nnz_basis_vec.begin(), nnz_basis_vec.end()
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
    b(i) = final_conditions.at(i - (num_rows - num_final)).second;

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

template <typename T, BoundaryCondition BC, Extrapolation EXT>
control_points::ControlPoints<T, BC> dense_nonuniform(
    size_t degree,
    knots::Knots<T, Curve::NON_UNIFORM, BC, EXT> const &knots,
    std::vector<T> const &x,
    std::vector<T> const &y,
    std::vector<std::pair<size_t, T>> const &initial_conditions,
    std::vector<std::pair<size_t, T>> const &final_conditions
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

  size_t num_initial  = initial_conditions.size();
  size_t num_final    = final_conditions.size();
  size_t num_rows     = x.size() + num_initial + num_final;
  Eigen::MatrixX<T> A = Eigen::MatrixX<T>::Zero(num_rows, num_cols);
  Eigen::VectorX<T> b;
  b.reserve(num_rows);

  size_t index{0}, i{0};
  std::vector<T> nnz_basis_vec(degree + 1, (T)0);
  T initial_x = x.front();
  T final_x   = x.back();
  for (; i < num_initial; i++)
  {
    // TODO: add degree of derivative, same for the others.
    // TODO: same about the knots.find.
    index = bspline::BSpline<T, Curve::NON_UNIFORM, BC, EXT>::nnz_basis(
        degree, knots, initial_x, nnz_basis_vec.begin(), nnz_basis_vec.end()
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
    b(i) = initial_conditions.at(i).second;

    std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
  }

  for (; i < x.size(); i++)
  {
    index = bspline::BSpline<T, Curve::NON_UNIFORM, BC, EXT>::nnz_basis(
        degree, knots, x.at(i - num_initial), nnz_basis_vec.begin(), nnz_basis_vec.end()
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
    b(i) = y.at(i - num_initial);

    std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
  }

  for (; i < num_rows; i++)
  {
    index = bspline::BSpline<T, Curve::NON_UNIFORM, BC, EXT>::nnz_basis(
        degree, knots, final_x, nnz_basis_vec.begin(), nnz_basis_vec.end()
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
    b(i) = final_conditions.at(i - (num_rows - num_final)).second;

    std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
  }

  Eigen::VectorX<T> res = A.colPivHouseholderQr().solve(b);

  return control_points::ControlPoints<T, BC>{
      {{res.data(), res.data() + res.rows() * res.cols()}}, degree
  };
}

template <typename T, BoundaryCondition BC, Extrapolation EXT>
control_points::ControlPoints<T, BC> sparse_nonuniform(
    size_t degree,
    knots::Knots<T, Curve::NON_UNIFORM, BC, EXT> const &knots,
    std::vector<T> const &x,
    std::vector<T> const &y,
    std::initializer_list<std::pair<size_t, T>> const &initial_conditions,
    std::initializer_list<std::pair<size_t, T>> const &final_conditions
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

  size_t num_initial = initial_conditions.size();
  size_t num_final   = final_conditions.size();
  size_t num_rows    = x.size() + num_initial + num_final;
  Eigen::SparseMatrix<T> A(num_rows, num_cols);
  A.reserve(num_cols * (degree + 1));
  Eigen::VectorX<T> b;
  b.reserve(num_rows);

  size_t index{0}, i{0};
  std::vector<T> nnz_basis_vec(degree + 1, (T)0);
  T initial_x = x.front();
  T final_x   = x.back();
  for (; i < num_initial; i++)
  {
    // TODO: add degree of derivative, same as before.
    // TODO: same about knots.find.
    index = bspline::BSpline<T, Curve::NON_UNIFORM, BC, EXT>::nnz_basis(
        degree, knots, initial_x, nnz_basis_vec.begin(), nnz_basis_vec.end()
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
    b(i) = initial_conditions.at(i).second;

    std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
  }

  for (; i < x.size(); i++)
  {
    index = bspline::BSpline<T, Curve::NON_UNIFORM, BC, EXT>::nnz_basis(
        degree, knots, x.at(i - num_initial), nnz_basis_vec.begin(), nnz_basis_vec.end()
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
    b(i) = y.at(i - num_initial);

    std::fill(nnz_basis_vec.begin(), nnz_basis_vec.end(), (T)0);
  }

  for (; i < num_rows; i++)
  {
    index = bspline::BSpline<T, Curve::NON_UNIFORM, BC, EXT>::nnz_basis(
        degree, knots, final_x, nnz_basis_vec.begin(), nnz_basis_vec.end()
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
    b(i) = final_conditions.at(i - (num_rows - num_final)).second;

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

} // namespace bsplinex::interpolate

#endif
