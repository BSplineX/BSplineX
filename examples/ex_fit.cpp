#include "BSplineX/bsplinex.hpp"
#include <iostream>

int main()
{
  // points to fit
  std::vector<double> const x_fit{2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9};
  std::vector<double> const y_fit{11.2, 22.3, 13.4, 14.5, 25.6, 36.7, 17.8};

  // initialize a cubic, uniform, periodic B-spline curve where knots are [0.1, 1.1,
  // ..., 10.1, 11.1]
  constexpr size_t degree{3};
  constexpr double knots_begin{0.1}, knots_end{12.0};
  constexpr size_t num_knots{12};

  auto bspline = bsplinex::factory::periodic_uniform(degree, knots_begin, knots_end, num_knots);

  // fit the curve to the points
  bspline.fit(x_fit, y_fit);

  // evaluate the curve at some points. Since the curve is periodic, the evaluation can done at any
  // point
  std::vector<double> const eval_x{-10.5, 3.0, 3.4, 5.1, 6.2, 15.0, 22.0};
  for (auto const &x : eval_x)
  {
    std::cout << "bspline.evaluate(" << x << ") = " << bspline.evaluate(x) << std::endl;
  }
}
