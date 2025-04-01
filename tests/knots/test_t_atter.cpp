// Standard includes
#include <cstddef>
#include <vector>

// Third-party includes
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// BSplineX includes
#include "BSplineX/knots/t_atter.hpp"
#include "BSplineX/knots/t_data.hpp"
#include "BSplineX/types.hpp"
#include "matchers.hpp"

using namespace Catch::Matchers;
using namespace bsplinex;
using namespace bsplinex::knots;

TEST_CASE("knots::Atter<T, C, BC> atter{knots::Data<T, C> data, degree}", "[t_atter]")
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  size_t const n = data_vec.size();
  Data<double, Curve::NON_UNIFORM> const data{data_vec};
  size_t constexpr degree{3};
  Atter<double, Curve::NON_UNIFORM, BoundaryCondition::PERIODIC> const atter{data, degree};

  SECTION("atter.size()") { REQUIRE(atter.size() == n + 2 * degree); }
  SECTION("atter.at(...)")
  {
    REQUIRE_THAT(
        atter.at(0), WithinAbsRel(data_vec.at(0) - (data_vec.at(n - 1) - data_vec.at(n - 4)))
    );
    REQUIRE_THAT(
        atter.at(1), WithinAbsRel(data_vec.at(0) - (data_vec.at(n - 1) - data_vec.at(n - 3)))
    );
    REQUIRE_THAT(
        atter.at(2), WithinAbsRel(data_vec.at(0) - (data_vec.at(n - 1) - data_vec.at(n - 2)))
    );
    for (size_t i{0}; i < n; i++)
    {
      REQUIRE(atter.at(i + degree) == data.at(i));
    }
    REQUIRE_THAT(
        atter.at(n + degree), WithinAbsRel(data_vec.at(n - 1) + (data_vec.at(1) - data_vec.at(0)))
    );
    REQUIRE_THAT(
        atter.at(n + degree + 1),
        WithinAbsRel(data_vec.at(n - 1) + (data_vec.at(2) - data_vec.at(0)))
    );
    REQUIRE_THAT(
        atter.at(n + degree + 2),
        WithinAbsRel(data_vec.at(n - 1) + (data_vec.at(3) - data_vec.at(0)))
    );
  }
  SECTION("atter.begin()")
  {
    auto it = atter.begin();
    REQUIRE(*it == atter.at(0));
  }
  SECTION("atter.end()")
  {
    auto it = atter.end();
    --it;
    REQUIRE(*it == atter.at(n + degree + 2));
  }
}

TEST_CASE("knots::Atter::iterator", "[t_atter]")
{
  std::vector<double> const data_vec{0.1, 1.3, 2.2, 4.9, 13.2};
  Data<double, Curve::NON_UNIFORM> const data{data_vec};
  size_t constexpr degree{3};
  Atter<double, Curve::NON_UNIFORM, BoundaryCondition::PERIODIC> const atter{data, degree};

  SECTION("*iterator")
  {
    auto it = atter.begin();
    REQUIRE(*it == atter.at(0));
  }
  SECTION("iterator++")
  {
    auto it = atter.begin();
    REQUIRE(*(it++) == atter.at(0));
    REQUIRE(*it == atter.at(1));
  }
  SECTION("++iterator")
  {
    auto it = atter.begin();
    REQUIRE(*(++it) == atter.at(1));
    REQUIRE(*it == atter.at(1));
  }
  SECTION("iterator--")
  {
    auto it = atter.end();
    --it;
    REQUIRE(*(it--) == atter.at(atter.size() - 1));
    REQUIRE(*it == atter.at(atter.size() - 2));
  }
  SECTION("--iterator")
  {
    auto it = atter.end();
    REQUIRE(*(--it) == atter.at(atter.size() - 1));
    REQUIRE(*it == atter.at(atter.size() - 1));
  }
  SECTION("iterator += 1")
  {
    auto it = atter.begin();
    REQUIRE(*(it += 1) == atter.at(1));
    REQUIRE(*it == atter.at(1));
  }
  SECTION("iterator + 1")
  {
    auto it = atter.begin();
    REQUIRE(*(it + 1) == atter.at(1));
    REQUIRE(*it == atter.at(0));
  }
  SECTION("iterator -= 1")
  {
    auto it = atter.end();
    REQUIRE(*(it -= 1) == atter.at(atter.size() - 1));
    REQUIRE(*it == atter.at(atter.size() - 1));
  }
  SECTION("iterator - 1")
  {
    auto it = atter.end();
    --it;
    REQUIRE(*(it - 1) == atter.at(atter.size() - 2));
    REQUIRE(*it == atter.at(atter.size() - 1));
  }
  SECTION("iterator - iterator")
  {
    auto it1 = atter.begin();
    auto it2 = atter.end();
    REQUIRE((it2 - it1) == (int)atter.size());
  }
  SECTION("iterator == iterator")
  {
    auto it1 = atter.begin();
    auto it2 = atter.begin();
    REQUIRE(it1 == it2);
  }
  SECTION("iterator != iterator")
  {
    auto it1 = atter.begin();
    auto it2 = atter.end();
    REQUIRE(it1 != it2);
  }
  SECTION("iterator[]")
  {
    auto it = atter.begin();
    REQUIRE(it[atter.size() - 2] == atter.at(atter.size() - 2));
  }
  SECTION("iterator > iterator")
  {
    auto it1 = atter.begin();
    auto it2 = atter.end();
    auto it3 = atter.begin();
    REQUIRE(it2 > it1);
    REQUIRE_FALSE(it1 > it2);
    REQUIRE_FALSE(it3 > it1);
  }
  SECTION("iterator < iterator")
  {
    auto it1 = atter.begin();
    auto it2 = atter.end();
    auto it3 = atter.begin();
    REQUIRE(it1 < it2);
    REQUIRE_FALSE(it2 < it1);
    REQUIRE_FALSE(it1 < it3);
  }
  SECTION("iterator >= iterator")
  {
    auto it1 = atter.begin();
    auto it2 = atter.end();
    auto it3 = atter.begin();
    REQUIRE(it2 >= it1);
    REQUIRE(it3 >= it1);
    REQUIRE_FALSE(it1 >= it2);
  }
  SECTION("iterator <= iterator")
  {
    auto it1 = atter.begin();
    auto it2 = atter.end();
    auto it3 = atter.begin();
    REQUIRE(it1 <= it2);
    REQUIRE(it1 <= it3);
    REQUIRE_FALSE(it2 <= it1);
  }
}
