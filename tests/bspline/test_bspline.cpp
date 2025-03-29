// Third-party includes
#include <algorithm>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

// BSplineX includes
#include "BSplineX/bspline/bspline_lsq.hpp"
#include "BSplineX/bspline/bspline_types.hpp"
#include "BSplineX/defines.hpp"
#include "BSplineX/types.hpp"

using namespace Catch::Matchers;
using namespace bsplinex;
using namespace bsplinex::bspline;
using namespace bsplinex::constants;

using real_t = double;

// clang-format off
#define BSPLINE_TEST_TYPES \
types::OpenUniform<real_t>, \
types::OpenUniformConstant<real_t>, \
types::OpenNonUniform<real_t>, \
types::OpenNonUniformConstant<real_t>, \
types::ClampedUniform<real_t>, \
types::ClampedUniformConstant<real_t>, \
types::ClampedNonUniform<real_t>, \
types::ClampedNonUniformConstant<real_t>, \
types::PeriodicUniform<real_t>, \
types::PeriodicNonUniform<real_t>

// clang-format on

static std::mt19937 &get_device()
{
  static std::random_device device{};
  static std::mt19937 rng{device()};

  rng.seed(05535);

  return rng;
}

static std::vector<real_t>
random_vector(std::mt19937 &rng, size_t num_elems, real_t min = -10.0, real_t max = 10.0)
{
  std::uniform_real_distribution unif{min, max};
  std::vector<real_t> vec(num_elems);
  std::generate(vec.begin(), vec.end(), [&unif, &rng]() { return unif(rng); });

  return vec;
}

static std::vector<real_t> uniform_vector(size_t num_elems, real_t min = -10.0, real_t max = 10.0)
{
  std::vector<real_t> vec(num_elems);
  real_t step_size = (max - min) / (num_elems - 1);
  std::generate(
      vec.begin(),
      std::prev(vec.end(), 1),
      [min, step_size, i = ZERO<size_t>]() mutable { return min + step_size * i++; }
  );

  // Avoid numerical errors
  vec.back() = max;

  return vec;
}

static std::vector<real_t>
sorted_random_vector(std::mt19937 &rng, size_t num_elems, real_t min = -10.0, real_t max = 10.0)
{
  std::vector<real_t> vec = random_vector(rng, num_elems, min, max);
  std::sort(vec.begin(), vec.end());

  return vec;
}

template <typename BSplineType>
BSplineType random_bspline(std::mt19937 &rng, size_t degree, size_t num_ctrl)
{
  std::vector<real_t> ctrl_pts = random_vector(rng, num_ctrl, -1.0, 1.0);

  size_t num_knots;
  if constexpr (BoundaryCondition::OPEN == BSplineType::boundary_condition_type)
  {
    num_knots = num_ctrl + degree + 1;
  }
  else if constexpr (BoundaryCondition::CLAMPED == BSplineType::boundary_condition_type)
  {
    num_knots = num_ctrl - degree + 1;
  }
  else if constexpr (BoundaryCondition::PERIODIC == BSplineType::boundary_condition_type)
  {
    num_knots = num_ctrl + 1;
  }
  else
  {
    throw std::runtime_error("Unknown BSpline Boundary Condition");
  }
  constexpr real_t knots_begin{-100.0};
  constexpr real_t knots_end{100.0};

  if constexpr (Curve::UNIFORM == BSplineType::curve_type)
  {
    BSplineType bspline{{knots_begin, knots_end, num_knots}, {ctrl_pts}, degree};

    return bspline;
  }
  else if constexpr (Curve::NON_UNIFORM == BSplineType::curve_type)
  {
    std::vector<real_t> knots = sorted_random_vector(rng, num_knots, knots_begin, knots_end);

    BSplineType bspline{{knots}, {ctrl_pts}, degree};

    return bspline;
  }
  else
  {
    throw std::runtime_error("Unkown BSpline Curve.");
  }
}

TEMPLATE_TEST_CASE("masinag", "[bspline][template][product]", BSPLINE_TEST_TYPES)
{
  auto &rng = get_device();

  size_t const degree =
      GENERATE(1); // TODO: test with bigger degrees once we implement derivatives for interpolation

  using ctrl_val_t              = std::pair<size_t, size_t>;
  auto const [ctrl_pts, values] = GENERATE(ctrl_val_t{50, 50}, ctrl_val_t{600, 2000});

  TestType bspline = random_bspline<TestType>(rng, degree, ctrl_pts);

  auto [x_min, x_max] = bspline.domain();

  std::vector<real_t> x;
  if constexpr (Curve::UNIFORM == TestType::curve_type)
  {
    x = uniform_vector(values, x_min, x_max);
  }
  else if constexpr (Curve::NON_UNIFORM == TestType::curve_type)
  {
    x = sorted_random_vector(rng, values, x_min, x_max);
  }
  else
  {
    throw std::runtime_error("Unkown BSpline Curve.");
  }

  std::vector<real_t> y = bspline.evaluate(x);

  SECTION("evaluate(...)")
  {
    for (size_t i{0}; i < x.size(); i++)
    {
      REQUIRE_THAT(bspline.evaluate(x.at(i)), WithinRel(y.at(i)));
    }
  }

  if constexpr (Extrapolation::CONSTANT == TestType::extrapolation_type)
  {
    SECTION("evaluate(...) - extrapolate constant")
    {
      for (size_t i{0}; i < 50; i++)
      {
        REQUIRE_THAT(
            bspline.evaluate(x_min - static_cast<real_t>(i)), WithinRel(bspline.evaluate(x_min))
        );
        REQUIRE_THAT(
            bspline.evaluate(x_max + static_cast<real_t>(i)), WithinRel(bspline.evaluate(x_max))
        );
      }
    }
  }
  else if constexpr (Extrapolation::PERIODIC == TestType::extrapolation_type)
  {
    SECTION("evaluate(...) - extrapolate periodic")
    {
      real_t const period = x_max - x_min;
      for (size_t i{0}; i < x.size(); i++)
      {
        REQUIRE_THAT(bspline.evaluate(x.at(i) + period), WithinRel(y.at(i), 1e-8));
        REQUIRE_THAT(bspline.evaluate(x.at(i) - period), WithinRel(y.at(i), 1e-8));
      }
    }
  }

  SECTION("compute_basis(...)")
  {
    auto c_data = bspline.get_control_points();

    for (size_t i{0}; i < x.size(); i++)
    {
      std::vector<real_t> basis = bspline.basis(x.at(i));
      REQUIRE(basis.size() == c_data.size());
      real_t res{0.0};
      for (size_t j{0}; j < basis.size(); j++)
      {
        res += basis.at(j) * c_data.at(j);
      }
      REQUIRE_THAT(res, WithinRel(y.at(i), 1e-10));
    }
  }

  SECTION("fit(...)")
  {
    bspline.fit(x, y);

    for (size_t i{0}; i < x.size(); i++)
    {
      REQUIRE_THAT(bspline.evaluate(x.at(i)), WithinRel(y.at(i), 1e-10));
    }
  }

  SECTION("interpolate(...)")
  {
    std::vector<lsq::Condition<real_t>> additional;

    if constexpr (BoundaryCondition::PERIODIC != TestType::boundary_condition_type)
    {
      for (size_t i{0}; i < degree - 1; i++)
      {
        // TODO: should be done with a derivative to avoid singular matrices
        additional.emplace_back(x.at(degree + i), y.at(degree + i), 0);
      }
    }
    else
    {
      y.back() = y.front();
    }

    bspline.interpolate(x, y, additional);
    size_t padding{0};
    if constexpr (BoundaryCondition::OPEN == TestType::boundary_condition_type)
    {
      padding = degree;
    }
    for (size_t i{padding}; i < x.size() - padding; i++)
    {
      REQUIRE_THAT(bspline.evaluate(x.at(i)), WithinRel(y.at(i), 1e-8));
    }
  }

  if constexpr (BoundaryCondition::PERIODIC != TestType::boundary_condition_type)
  {
    SECTION("interpolate(...) - invalid additional conditions")
    {
      std::vector<lsq::Condition<real_t>> additional;
      size_t const additional_size = degree - 1;
      for (size_t i{0}; i < additional_size; i++)
      {
        additional.emplace_back(x.back() + ONE<real_t>, y.back(), 0);
      }
      if (additional_size > 0)
      {
        REQUIRE_THROWS_WITH(
            bspline.interpolate(x, y, additional),
            "Additional conditions must lie inside the knots interval."
        );
      }
    }
  }
}
