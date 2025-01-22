#include "BSplineX/bsplinex.hpp"
#include <cmath>
#include <iostream>
#include <vector>

int main()
{
  // points to interpolate
  const std::vector<double> x_interp{2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 40.2};
  std::vector<double> y_interp{11.2, 22.3, 13.4, 14.5, 25.6, 36.7, 17.8, 11.2};

  constexpr std::size_t degree{3};
  const std::vector<double> knots{0.1, 1.3, 2.2, 2.2, 4.9, 6.3, 6.3, 6.3, 13.2};

  auto bspline = bsplinex::factory::periodic_nonuniform(degree, knots);

  // interpolate the curve to the points
  bspline.interpolate(x_interp, y_interp, {});

  double x, y;
  for (size_t i{0}; i < x_interp.size(); i++)
  {
    x = x_interp.at(i);
    y = y_interp.at(i);
    std::cout << "bspline.evaluate(" << x << ") = " << bspline.evaluate(x) << " == " << y
              << " -> abs(error) = " << std::abs(bspline.evaluate(x) - y) << std::endl;
  }
}
