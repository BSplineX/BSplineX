// Third-party includes
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <random>

// BSplineX includes
#include "BSplineX/bspline/bspline_types.hpp"

using namespace Catch::Matchers;
using namespace bsplinex;
using namespace bsplinex::bspline;

TEST_CASE(
    "bspline::BSpline<T, C, BoundaryCondition::OPEN, Extrapolation::NONE> "
    "bspline{knots::Data<T, C> t_data, "
    "control_points::Data<T> c_data, degree}",
    "[bspline]"
)
{
  size_t degree{3};

  std::vector<double> t_data_vec{0.1, 1.3, 2.2, 2.2, 4.9, 6.3, 6.3, 6.3, 13.2};
  knots::Data<double, Curve::NON_UNIFORM> t_data{t_data_vec};

  std::vector<double> c_data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  control_points::Data<double> c_data{c_data_vec};

  types::OpenNonUniform<double> bspline{t_data, c_data, degree};

  // clang-format off
  std::vector<double> x_values{2.2, 2.3000000000000003, 2.4000000000000004, 2.5000000000000004, 2.6000000000000005, 2.7000000000000006, 2.8000000000000007, 2.900000000000001, 3.000000000000001, 3.100000000000001, 3.200000000000001, 3.300000000000001, 3.4000000000000012, 3.5000000000000013, 3.6000000000000014, 3.7000000000000015, 3.8000000000000016, 3.9000000000000017, 4.000000000000002, 4.100000000000001, 4.200000000000002, 4.3000000000000025, 4.400000000000002, 4.500000000000002, 4.600000000000002, 4.700000000000003, 4.8000000000000025, 4.900000000000002, 5.000000000000003, 5.100000000000003, 5.200000000000003, 5.3000000000000025, 5.400000000000003, 5.5000000000000036, 5.600000000000003, 5.700000000000003, 5.800000000000003, 5.900000000000004, 6.0000000000000036, 6.100000000000003, 6.200000000000004};
  std::vector<double> y_values{0.4, 0.4987905929445725, 0.5953834608104188, 0.6901102371457323, 0.7833025554987059, 0.8752920494175335, 0.9664103524504085, 1.0569890981455239, 1.1473599200510731, 1.2378544517152497, 1.3288043266862466, 1.420541178512258, 1.5133966407414765, 1.6077023469220952, 1.7037899306023085, 1.8019910253303089, 1.9026372646542904, 2.0060602821224456, 2.112591711282968, 2.2225631856840518, 2.3363063388738903, 2.4541528044006755, 2.5764342158126015, 2.7034822066578617, 2.83562841048465, 2.9732044608411603, 3.1165419912755823, 3.2659726353361123, 3.4243850625148546, 3.604896086079547, 3.8231795552418366, 4.094909319213373, 4.435759227205808, 4.861403128430789, 5.387514872099959, 6.029768307424969, 6.8038372836174785, 7.725395649889126, 8.810117255451555, 10.073675949516419, 11.53174558129538};
  // clang-format on

  SECTION("bspline.evaluate(...)")
  {
    for (size_t i{0}; i < x_values.size(); i++)
    {
      REQUIRE_THAT(bspline.evaluate(x_values.at(i)), WithinRel(y_values.at(i)));
    }
  }

  SECTION("bspline.compute_basis(...)")
  {
    for (size_t i{0}; i < x_values.size(); i++)
    {
      std::vector<double> basis = bspline.basis(x_values.at(i));
      REQUIRE(basis.size() == c_data.size());
      double res{0.0};
      for (size_t j{0}; j < basis.size(); j++)
      {
        res += basis.at(j) * c_data.at(j);
      }
      REQUIRE_THAT(res, WithinRel(y_values.at(i)));
    }
  }

  SECTION("bspline.fit(...) dense")
  {
    bspline.fit(x_values, y_values);
    auto control_points = bspline.get_control_points();
    for (size_t i{0}; i < c_data.size(); i++)
    {
      REQUIRE_THAT(control_points.at(i), WithinRel(c_data.at(i), 1e-6));
    }
    for (size_t i{0}; i < x_values.size(); i++)
    {
      REQUIRE_THAT(bspline.evaluate(x_values.at(i)), WithinRel(y_values.at(i), 1e-6));
    }
  }

  SECTION("bspline.fit(...) sparse")
  {
    // Prepare a normal distribution
    std::random_device device{};
    std::mt19937 rng{device()};
    rng.seed(05535);
    std::normal_distribution norm{0.0, 1.0};

    // Generated big knots and ctrl points
    std::vector<double> big_ctrl_pts(513);
    std::vector<double> big_knots(big_ctrl_pts.size() + 3 + 1);
    std::generate(big_ctrl_pts.begin(), big_ctrl_pts.end(), [&norm, &rng]() { return norm(rng); });
    std::generate(big_knots.begin(), big_knots.end(), [n = 0]() mutable { return (double)n++; });

    types::OpenNonUniform<double> big_bspline{big_knots, big_ctrl_pts, degree};

    // Prepare a uniform distribution
    std::uniform_real_distribution unif{big_knots.at(3), big_knots.at(big_knots.size() - 4)};

    // Randomly sample points
    std::vector<double> big_x(2000);
    std::vector<double> big_y(big_x.size());
    std::generate(big_x.begin(), big_x.end(), [&unif, &rng]() { return unif(rng); });
    std::generate(
        big_y.begin(),
        big_y.end(),
        [i = 0, &big_bspline, &big_x]() mutable { return big_bspline.evaluate(big_x.at(i++)); }
    );

    big_bspline.fit(big_x, big_y);
    auto control_points = big_bspline.get_control_points();
    for (size_t i{0}; i < big_ctrl_pts.size(); i++)
    {
      REQUIRE_THAT(control_points.at(i), WithinRel(big_ctrl_pts.at(i), 1e-6));
    }
    for (size_t i{0}; i < big_x.size(); i++)
    {
      REQUIRE_THAT(big_bspline.evaluate(big_x.at(i)), WithinRel(big_y.at(i), 1e-6));
    }
  }
}

TEST_CASE(
    "bspline::BSpline<T, C, BoundaryCondition::OPEN, Extrapolation::CONSTANT> "
    "bspline{knots::Data<T, C> t_data, "
    "control_points::Data<T> c_data, degree}",
    "[bspline]"
)
{
  size_t degree{3};

  std::vector<double> t_data_vec{0.1, 1.3, 2.2, 2.2, 4.9, 6.3, 6.3, 6.3, 13.2};
  knots::Data<double, Curve::NON_UNIFORM> t_data{t_data_vec};

  std::vector<double> c_data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  control_points::Data<double> c_data{c_data_vec};

  types::OpenNonUniformConstant<double> bspline{t_data, c_data, degree};

  // clang-format off
  std::vector<double> x_values{-1.0, 0.0, 1.0, 2.2, 3.225, 4.25, 5.275, 6.3, 6.4, 6.9, 7.1};
  std::vector<double> y_values{0.4, 0.4, 0.4, 0.4, 1.3516518061271148, 2.3946959304983997, 4.0211091244533534, 13.2, 13.2, 13.2, 13.2};
  // clang-format on

  SECTION("bspline.evaluate(...)")
  {
    for (size_t i{0}; i < x_values.size(); i++)
    {
      REQUIRE_THAT(bspline.evaluate(x_values.at(i)), WithinRel(y_values.at(i)));
    }
  }
}
