#include <fstream>

// Third-party includes
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

// BSplineX includes

using namespace Catch::Matchers;

TEST_CASE("test", "[test]")
{
  std::ifstream f("example.json");
  json data = json::parse(f);
}
