#include "BSplineX/bsplinex.hpp"
#include <iostream>

int main()
{
  // initialize a cubic, non-uniform, open B-spline curve with given knots and control points
  constexpr std::size_t degree{3};
  const std::vector<double> knots{0.1, 1.3, 2.2, 2.2, 4.9, 6.3, 6.3, 6.3, 13.2};
  const std::vector<double> control_points{0.1, 1.3, 2.2, 4.9, 13.2};

  auto bspline = bsplinex::factory::open_nonuniform(degree, knots, control_points);

  // evaluate the curve at some points
  std::vector<double> eval_x{3.0, 3.4, 5.1, 6.2};
  for (const double x : eval_x)
  {
    std::cout << "bspline.evaluate(" << x << ") = " << bspline.evaluate(x) << std::endl;
  }
}
