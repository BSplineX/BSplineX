// Standard includes
#include <cmath>
#include <string>

// Third-party includes
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/catch_test_macros.hpp>

// BSplineX includes
#include "BSplineX/bspline/bspline.hpp"

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
  knots::Data<double, Curve::UNIFORM> t_data{0.1, 13.2, static_cast<size_t>(9)};
  control_points::Data<double> c_data{{0.1, 1.3, 2.2, 4.9, 13.2}};

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
    double step = (stop - start) / (double)steps;
    vec.resize(steps);
    for (size_t i{0}; i < steps; i++)
    {
      vec.at(i) = start + (double)i * step;
    }
  };

  size_t knots_num{0};
  std::vector<double> ctrl_pts{};
  std::vector<double> x_data{};
  std::vector<double> y_data{};
  double res{0.0};
  double eval_elems{0.0};
  double start{45.0};
  double stop{49.0};
  for (size_t j{3}; j < 12; j++)
  {
    knots_num = (size_t)std::pow(2.0, j);
    fill(start, stop, knots_num - degree - 1, ctrl_pts);
    BSpline<double, Curve::UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE> bspline{
        {0.0, 100.0, knots_num}, {ctrl_pts}, degree
    };

    for (size_t i{4}; i < 5; i++)
    {
      eval_elems = std::pow(10.0, i);
      fill(start, stop, static_cast<size_t>(eval_elems), x_data);
      res = 0.0;

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
      for (auto x : x_data)
      {
        y_data.push_back(bspline.evaluate(x));
      }

      BENCHMARK(
          "bspline.fit - knots: " + std::to_string(knots_num) +
          " points: " + std::to_string((size_t)eval_elems)
      )
      {
        return bspline.fit(x_data, y_data);
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
  knots::Data<double, Curve::NON_UNIFORM> t_data{{0.1, 1.3, 2.2, 2.2, 4.9, 6.3, 6.3, 6.3, 13.2}};
  control_points::Data<double> c_data{{0.1, 1.3, 2.2, 4.9, 13.2}};

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

  BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE> bspline{
      t_data, c_data, degree
  };

  auto fill = [](double start, double stop, size_t steps, std::vector<double> &vec)
  {
    double step = (stop - start) / (double)steps;
    vec.resize(steps);
    for (size_t i{0}; i < steps; i++)
    {
      vec.at(i) = start + (double)i * step;
    }
  };

  size_t knots_num{0};
  std::vector<double> knots{};
  std::vector<double> ctrl_pts{};
  std::vector<double> x_data{};
  double res{0.0};
  double eval_elems{0.0};
  double start{45.0};
  double stop{49.0};
  for (size_t j{3}; j < 11; j++)
  {
    knots_num = (size_t)std::pow(2.0, j);
    fill(0.0, 100.0, knots_num, knots);
    fill(start, stop, knots_num - degree - 1, ctrl_pts);
    BSpline<double, Curve::NON_UNIFORM, BoundaryCondition::OPEN, Extrapolation::NONE> bspline{
        {knots}, {ctrl_pts}, degree
    };

    for (size_t i{4}; i < 5; i++)
    {
      eval_elems = std::pow(10.0, i);
      fill(start, stop, static_cast<size_t>(eval_elems), x_data);
      res = 0.0;
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
