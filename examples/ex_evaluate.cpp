#include "BSplineX/bsplinex.hpp"
#include <cstddef>
#include <iostream>
#include <vector>

int main()
{
  // initialize a cubic, non-uniform, open B-spline curve with given knots and control points
  constexpr size_t degree{3};
  std::vector<double> const knots{0.1, 1.3, 2.2, 2.2, 4.9, 6.3, 6.3, 6.3, 13.2};
  std::vector<double> const control_points{0.1, 1.3, 2.2, 4.9, 13.2};

  auto bspline = bsplinex::open_nonuniform(degree, knots, control_points);

  // evaluate the curve at some points
  std::vector<double> const eval_x{3.0, 3.4, 5.1, 6.2};
  for (double const x : eval_x)
  {
    std::cout << "bspline.evaluate(" << x << ") = " << bspline.evaluate(x) << "\n";
  }
}
