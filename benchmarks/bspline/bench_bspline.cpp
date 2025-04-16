// Standard includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <string>
#include <vector>

// Third-party includes
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/benchmark/catch_chronometer.hpp>
#include <catch2/benchmark/catch_constructor.hpp>
#include <catch2/catch_test_macros.hpp>

// BSplineX includes
#include "BSplineX/bspline/bspline.hpp"
#include "BSplineX/control_points/c_data.hpp"
#include "BSplineX/knots/t_data.hpp"
#include "BSplineX/types.hpp"

using namespace bsplinex;
using namespace bsplinex::bspline;

/**
 * Rerefence data can be generated using the `reference_sbpline.py` file in the
 * tests folder
 */

TEST_CASE(
    "benchmark bspline::BSpline<double, Curve::UNIFORM, "
    "BoundaryCondition::OPEN, Extrapolation::NONE>",
    "[bspline]"
)
{
  size_t degree{3};
  knots::Data<double, Curve::UNIFORM> const t_data{0.1, 13.2, static_cast<size_t>(9)};
  control_points::Data<double> const c_data{{0.1, 1.3, 2.2, 4.9, 13.2}};

  BENCHMARK_ADVANCED("bspline.construct")(Catch::Benchmark::Chronometer meter)
  {
    std::vector<Catch::Benchmark::storage_for<
        BSpline<double, Curve::UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE>>>
        storage(meter.runs());
    meter.measure([&](int i) { storage[i].construct(t_data, c_data, degree); });
  };

  BENCHMARK_ADVANCED("bspline.destroy")(Catch::Benchmark::Chronometer meter)
  {
    std::vector<Catch::Benchmark::destructable_object<
        BSpline<double, Curve::UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE>>>
        storage(meter.runs());
    for (auto &&o : storage)
    {
      o.construct(t_data, c_data, degree);
    }

    meter.measure([&](int i) { storage[i].destruct(); });
  };

  auto fill = [](double start, double stop, size_t steps, std::vector<double> &vec)
  {
    double const step = (stop - start) / (double)steps;
    vec.resize(steps);
    for (size_t i{0}; i < steps; i++)
    {
      vec.at(i) = start + (double)i * step;
    }
  };

  std::vector<double> ctrl_pts{};
  std::vector<double> x_data{};
  std::vector<double> y_data{};
  double const start{45.0};
  double const stop{49.0};
  size_t const min_knots_pow{3};
  size_t const max_knots_pow{12};
  for (size_t j{min_knots_pow}; j < max_knots_pow; j++)
  {
    auto const knots_num = static_cast<size_t>(std::pow(2.0, j));
    fill(start, stop, knots_num - degree - 1, ctrl_pts);
    BSpline<double, Curve::UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE> bspline{
        knots::Data<double, Curve::UNIFORM>{0.0, 100.0, knots_num},
        control_points::Data<double>{ctrl_pts},
        degree
    };

    size_t const min_eval_pow{4};
    size_t const max_eval_pow{5};
    double const eval_pow_base{10.0};
    for (size_t i{min_eval_pow}; i < max_eval_pow; i++)
    {
      double const eval_elems{std::pow(eval_pow_base, i)};
      fill(start, stop, static_cast<size_t>(eval_elems), x_data);
      double res{0.0};

      BENCHMARK(
          "bspline.evaluate - knots: " + std::to_string(knots_num) +
          " evals: " + std::to_string(static_cast<size_t>(eval_elems))
      )
      {
        for (auto x : x_data)
        {
          res = bspline.evaluate(x);
        }
        return res;
      };

      y_data.clear();
      y_data.reserve(x_data.size());
      std::transform(
          x_data.begin(),
          x_data.end(),
          std::back_inserter(y_data),
          [&bspline](double x) { return bspline.evaluate(x); }
      );

      BENCHMARK(
          "bspline.fit - knots: " + std::to_string(knots_num) +
          " points: " + std::to_string((size_t)eval_elems)
      )
      {
        bspline.fit(x_data, y_data);
      };
    }
  }
}

TEST_CASE(
    "benchmark bspline::BSpline<double, Curve::NON_UNIFORM, "
    "BoundaryCondition::OPEN, Extrapolation::NONE>",
    "[bspline]"
)
{
  size_t degree{3};
  knots::Data<double, Curve::NON_UNIFORM> const t_data{
      {0.1, 1.3, 2.2, 2.2, 4.9, 6.3, 6.3, 6.3, 13.2}
  };
  control_points::Data<double> const c_data{{0.1, 1.3, 2.2, 4.9, 13.2}};

  BENCHMARK_ADVANCED("bspline.construct")(Catch::Benchmark::Chronometer meter)
  {
    std::vector<Catch::Benchmark::storage_for<
        BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE>>>
        storage(meter.runs());
    meter.measure([&](int i) { storage[i].construct(t_data, c_data, degree); });
  };

  BENCHMARK_ADVANCED("bspline.destroy")(Catch::Benchmark::Chronometer meter)
  {
    std::vector<Catch::Benchmark::destructable_object<
        BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE>>>
        storage(meter.runs());
    for (auto &&o : storage)
    {
      o.construct(t_data, c_data, degree);
    }

    meter.measure([&](int i) { storage[i].destruct(); });
  };

  auto fill = [](double start, double stop, size_t steps, std::vector<double> &vec)
  {
    double const step = (stop - start) / (double)steps;
    vec.resize(steps);
    for (size_t i{0}; i < steps; i++)
    {
      vec.at(i) = start + (double)i * step;
    }
  };

  std::vector<double> knots{};
  std::vector<double> ctrl_pts{};
  std::vector<double> x_data{};
  double const start{45.0};
  double const stop{49.0};
  size_t const min_knots_pow{3};
  size_t const max_knots_pow{11};
  for (size_t j{min_knots_pow}; j < max_knots_pow; j++)
  {
    auto const knots_num = static_cast<size_t>(std::pow(2.0, j));
    fill(0.0, 100.0, knots_num, knots);
    fill(start, stop, knots_num - degree - 1, ctrl_pts);
    BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE> bspline{
        knots::Data<double, Curve::NON_UNIFORM>{knots},
        control_points::Data<double>{ctrl_pts},
        degree
    };

    size_t const min_eval_pow{4};
    size_t const max_eval_pow{5};
    double const eval_pow_base{10.0};
    for (size_t i{min_eval_pow}; i < max_eval_pow; i++)
    {
      double const eval_elems{std::pow(eval_pow_base, i)};
      fill(start, stop, static_cast<size_t>(eval_elems), x_data);
      double res{0.0};
      BENCHMARK(
          "bspline.evaluate - knots: " + std::to_string(knots_num) +
          " evals: " + std::to_string(static_cast<size_t>(eval_elems))
      )
      {
        for (auto x : x_data)
        {
          res = bspline.evaluate(x);
        }
        return res;
      };
    }
  }
}
