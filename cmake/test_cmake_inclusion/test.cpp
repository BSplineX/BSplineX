#include "BSplineX/bsplinex.hpp"
#include <iostream>

int main()
{
  auto bspline = bsplinex::make_periodic_nonuniform(3, {0, 1, 2, 3, 4});

  constexpr double x{3.0};
  std::cout << bspline.evaluate(x) << "\n";

  return 0;
}
