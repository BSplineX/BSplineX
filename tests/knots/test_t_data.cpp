// Standard includes
#include <cstddef>
#include <stdexcept>
#include <vector>

// Third-party includes
#include <catch2/catch_test_macros.hpp>

// BSplineX includes
#include "BSplineX/knots/t_data.hpp"
#include "BSplineX/types.hpp"

using namespace bsplinex;
using namespace bsplinex::knots;

TEST_CASE("knots::Data<double, Curve::UNIFORM> data{begin, end, num_elems}", "[t_data]")
{
  double constexpr begin{0.0};
  double constexpr end{10.0};
  double constexpr step{2.5};
  size_t constexpr num_elems{5};
  Data<double, Curve::UNIFORM> const data{begin, end, num_elems};

  SECTION("data.size()") { REQUIRE(data.size() == num_elems); }
  SECTION("data.at(...)")
  {
    for (size_t i{0}; i < data.size(); i++)
    {
      REQUIRE(data.at(i) == i * step);
    }
  }
  SECTION("data.slice(...)")
  {
    std::vector<double> const slice = data.slice(1, 4);
    REQUIRE(slice.size() == 3);
    for (size_t i{0}; i < slice.size(); i++)
    {
      REQUIRE(slice.at(i) == data.at(i + 1));
    }
  }
}

TEST_CASE("knots::Data<double, Curve::UNIFORM> data{data_vec}", "[t_data]")
{
  SECTION("Non-uniform vector")
  {
    std::vector<double> const data_vec{0.0, 2.5, 5.0, 7.5, 9.0};
    auto test = [&data_vec]() { Data<double, Curve::UNIFORM> const _(data_vec); };
    REQUIRE_THROWS_AS(test(), std::runtime_error);
  }
  SECTION("Uniform descending vector")
  {
    std::vector<double> const data_vec{10.0, 7.5, 5.0, 2.5, 0.0};
    auto test = [&data_vec]() { Data<double, Curve::UNIFORM> const _(data_vec); };
    REQUIRE_THROWS_AS(test(), std::runtime_error);
  }
  std::vector<double> const data_vec{0.0, 2.5, 5.0, 7.5, 10.0};
  Data<double, Curve::UNIFORM> const data{data_vec};

  SECTION("data.size()") { REQUIRE(data.size() == data_vec.size()); }
  SECTION("data.at(...)")
  {
    for (size_t i{0}; i < data.size(); i++)
    {
      REQUIRE(data.at(i) == data_vec.at(i));
    }
  }
  SECTION("data.slice(...)")
  {
    std::vector<double> const slice = data.slice(1, 4);
    REQUIRE(slice.size() == 3);
    for (size_t i{0}; i < slice.size(); i++)
    {
      REQUIRE(slice.at(i) == data.at(i + 1));
    }
  }
}

TEST_CASE("knots::Data<double, Curve::NON_UNIFORM> data{data_vec}", "[t_data]")
{
  std::vector<double> const data_vec{0.0, 1.3, 2.2, 4.9, 13.2};
  Data<double, Curve::NON_UNIFORM> const data{data_vec};

  SECTION("data.size()") { REQUIRE(data.size() == data_vec.size()); }
  SECTION("data.at(...)")
  {
    for (size_t i{0}; i < data.size(); i++)
    {
      REQUIRE(data.at(i) == data_vec.at(i));
    }
  }
  SECTION("data.slice(...)")
  {
    std::vector<double> const slice = data.slice(1, 4);
    REQUIRE(slice.size() == 3);
    for (size_t i{0}; i < slice.size(); i++)
    {
      REQUIRE(slice.at(i) == data.at(i + 1));
    }
  }
}
