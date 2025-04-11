#include "BSplineX/bsplinex.hpp"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

int main()
{
  // points to interpolate
  std::vector<double> const x_interp{2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 40.2};
  std::vector<double> const y_interp{11.2, 22.3, 13.4, 14.5, 25.6, 36.7, 17.8, 11.2};
  // additional condition: the value of the cond_d-th order derivative at cond_x must be cond_y
  double const cond_x = 6.7;
  double const cond_y = 25.6;
  size_t const cond_d = 1;

  constexpr size_t degree{2};

  auto bspline = bsplinex::factory::open_nonuniform(degree);

  // interpolate the curve to the points
  bspline.interpolate(x_interp, y_interp, {{cond_x, cond_y, cond_d}});

  for (size_t i{degree}; i < x_interp.size() - degree; i++)
  {
    double const x = x_interp.at(i);
    double const y = y_interp.at(i);
    std::cout << "bspline.evaluate(" << x << ") = " << bspline.evaluate(x) << " == " << y
              << " -> abs(error) = " << std::abs(bspline.evaluate(x) - y) << "\n";
  }
  std::cout << "bspline.evaluate(" << cond_x << ", derivative_order=" << cond_d
            << ") = " << bspline.evaluate(cond_x, cond_d) << " == " << cond_y
            << " -> abs(error) = " << std::abs(bspline.evaluate(cond_x, cond_d) - cond_y) << "\n";
}
