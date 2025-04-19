# BSplineX

The B-Spline library you have always dreamed of where performance and easy-to-use interfaces live in harmony.

## Features

This is a non-exhaustive list of the currently supported features.

- 1D B-Splines
  - Uniform/Non-uniform knots
  - Open/Clamped/Periodic boundary conditions
  - None/Constant/Periodic extrapolation
  - Least-squares fitting of the control points

## Installation

### Quick and dirty

`BSplineX` is a header-only library that depends only on `Eigen`, so the quick and dirty way of installing it is by simply copying the `include` directory to your project and make sure to have `Eigen` available however you see fit. Alternatively, you can do things properly and use `CMake`.

### CMake

If you are using CMake, you can add the library as a subdirectory in your project.

```cmake
add_subdirectory(path/to/BSplineX)
target_link_libraries(your_target PRIVATE BSplineX::BSplineX)
```

You can use `FetchContent` to download the library from GitHub.

```cmake
include(FetchContent)

# Optionally specify a custom path to fetch content to
set(FETCHCONTENT_BASE_DIR "path/to/your/dependencies")
fetchcontent_declare(
  BSplineX
  GIT_REPOSITORY https://github.com/BSplineX/BSplineX.git
  GIT_TAG        main
)
fetchcontent_makeavailable(BSplineX)
target_link_libraries(your_target PRIVATE BSplineX::BSplineX)
```

If you already have `BSplineX` somewhere on your system, you can use `find_pacakge` directly.

```cmake
# Optionally specify a custom path to find content from
list(APPEND CMAKE_PREFIX_PATH "path/to/your/dependencies")
find_package(
  BSplineX
  ${YOUR_DESIRED_BSPLINEX_VERSION}
  NO_MODULE
)

target_link_libraries(your_target PRIVATE BSplineX::BSplineX)
```

Since we are nice people, we also show you how to conditionally use `FetchContent` based if you already have the library or not.

```cmake
# Optionally specify a custom path to find content from
list(APPEND CMAKE_PREFIX_PATH "path/to/your/dependencies")
find_package(
  BSplineX
  ${YOUR_DESIRED_BSPLINEX_VERSION}
  NO_MODULE
)

if(NOT TARGET BSplineX::BSplineX)
  include(FetchContent)

  # Optionally specify a custom path to fetch content to
  set(FETCHCONTENT_BASE_DIR "path/to/your/dependencies")
  fetchcontent_declare(
      BSplineX
      GIT_REPOSITORY https://github.com/BSplineX/BSplineX.git
      GIT_TAG        main
  )

  fetchcontent_makeavailable(BSplineX)
endif()

target_link_libraries(your_target PRIVATE BSplineX::BSplineX)
```

## Usage

`BSplineX` is extremely easy to use, everything that is complex about B-Splines is hidden behind the curtain. Hereafter some quick examples.

### Creating a 1D B-Spline

```cpp
#include <iostream>
#include <vector>
#include "BSplineX/bsplinex.hpp"

int main()
{
  // Create a uniform clamped B-Spline with a 10 knots sequence in [0.0, 10.0]
  // with knots.size() - degree - 1 control points
  constexpr size_t degree{3};
  constexpr size_t num_knots{10};
  constexpr double min_knot{0.0};
  constexpr double max_knot{10.0};
  std::vector<double> ctrl_points{1.0, 0.3, 0.6, 0.9, 0.28, 9.32};

  auto bspline = bsplinex::clamped_uniform_constant(
    degree, min_knot, max_knot, num_knots, ctrl_points
  );

  // Now you can simply evaluate the B-Spline
  double x{3.26};
  std::cout << "B-Spline value at " << x << " is " << bspline.evaluate(x) << std::endl;

  return 0;
}
```

If you want to make it periodic or open, just change the template parameter, easy. Want to pass in a non-uniform knot sequence? Piece of cake, just set the `Curve` template parameter to `Curve::NON_UNIFORM` and pass in your knots vector. Boom, done, simple.

Now I know what you are thinking: "Dude, cool, but how on earth do I know which control points to use?". Fear not, and follow along ;).

### Fitting a 1D B-Spline

In the previous example, the knots are easy to understand: they are you discretisation grid, so you just have to choose how you want to discretise along that dimension. Control points on the other hand are not easily graspable. Fortunately you don't really need to know them a priori, you can just interpolate (coming soon!) or fit some data points. Let's see how, assuming we start from the previous example.

```cpp
...

  // Create the B-Spline without specifying the control points
  auto bspline = bsplinex::clamped_uniform_constant(
    degree, min_knot, max_knot, num_knots
  );

  // Assume somewhere you have loaded x- and y-values to fit you B-Spline to

  // Fit the B-Spline
  bspline.fit(x, y);

  // Now the B-Spline magically contains the best control points
  // (in a least-squares sense) that fit your data

...
```

### More examples

See the [examples](examples) folder for some more examples on how to use the library.

## Performance

`BSplineX` is all about performance. In fact, aside from being a cool and useful project, it is our playground to understand modern code optimisation techniques. Our goal is to be as fast as humanly possible (albeit without having to re-write the library by hand in assembly). Hereafter a non-exhaustive quick table with some performance results. If you want to can run the benchmarks on your own PC.
