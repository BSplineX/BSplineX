#include <array>
#include <numeric>
#include <type_traits>
#include <vector>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>

#include "BSplineX/defines.hpp"
#include "BSplineX/views.hpp"

using test_t = int;
constexpr size_t NUM_ELEMS{20};

// clang-format off
#define VIEW_TEST_TYPES \
decltype(std::vector<test_t>{}), \
decltype(std::array<test_t, NUM_ELEMS>{})
// clang-format on

template <typename T>
struct is_vector : std::false_type
{
};

template <typename T, typename A>
struct is_vector<std::vector<T, A>> : std::true_type
{
};

template <typename T>
inline constexpr bool is_vector_v = is_vector<T>::value;

using namespace bsplinex::views;

TEMPLATE_TEST_CASE("ArrayView", "[ArrayView][template][product]", VIEW_TEST_TYPES)
{
  using vec_type   = TestType;
  using iter       = typename vec_type::iterator;
  using const_iter = typename vec_type::const_iterator;

  vec_type data;
  if constexpr (is_vector_v<vec_type>)
  {
    data.resize(NUM_ELEMS);
  }
  std::iota(data.begin(), data.end(), bsplinex::constants::ONE<test_t>);
  vec_type const const_data = data;

  REQUIRE(data.size() == NUM_ELEMS);

  ArrayView<iter> view(data.begin(), data.end());
  ArrayView<const_iter> const_view(const_data.begin(), const_data.end());

  SECTION("Empty constructor") { [[maybe_unused]] ArrayView<iter> _view{}; }
  SECTION("Constructor") { ArrayView<iter> _view(data.begin(), data.end()); }

  SECTION("Size") { REQUIRE(view.size() == data.size()); }

  SECTION("Element access with operator[]")
  {
    for (size_t i{0}; i < view.size(); i++)
    {
      REQUIRE(view[i] == data[i]);
    }

    constexpr test_t new_elem{static_cast<test_t>(10)};
    view[1] = new_elem;
    REQUIRE(data[1] == new_elem);
  }

  SECTION("Element access with at()")
  {
    for (size_t i{0}; i < view.size(); i++)
    {
      REQUIRE(view.at(i) == data.at(i));
    }

    constexpr test_t new_elem{static_cast<test_t>(20)};
    view.at(3) = new_elem;
    REQUIRE(data.at(3) == new_elem);
  }

  SECTION("Front and back methods")
  {
    REQUIRE(view.front() == data.front());
    REQUIRE(view.back() == data.back());

    constexpr test_t new_front{static_cast<test_t>(100)};
    constexpr test_t new_back{static_cast<test_t>(500)};
    view.front() = new_front;
    view.back()  = new_back;
    REQUIRE(data.front() == new_front);
    REQUIRE(data.back() == new_back);
  }

  SECTION("Iterators")
  {
    size_t index{0};
    for (auto const &elem : view)
    {
      REQUIRE(elem == data.at(index++));
    }
  }

  SECTION("Const element access")
  {
    REQUIRE(const_view[0] == const_data[0]);
    REQUIRE(const_view.at(2) == const_data.at(2));
    REQUIRE(const_view.front() == const_data.front());
    REQUIRE(const_view.back() == const_data.back());
  }

  SECTION("Copy constructor")
  {
    ArrayView<iter> new_view{view};
    REQUIRE(new_view.size() == view.size());
    REQUIRE(new_view.front() == view.front());

    constexpr test_t new_elem{static_cast<test_t>(100)};
    new_view.front() = new_elem;
    REQUIRE(view.front() == new_elem);
    REQUIRE(data.front() == new_elem);
  }

  SECTION("Copy assignment operator")
  {
    ArrayView<iter> new_view{};
    new_view = view;
    REQUIRE(new_view.size() == view.size());
    REQUIRE(new_view.front() == view.front());

    constexpr test_t new_elem{static_cast<test_t>(200)};
    new_view.front() = new_elem;
    REQUIRE(view.front() == new_elem);
    REQUIRE(data.front() == new_elem);
  }

  SECTION("Move constructor")
  {
    auto tmp_view = view;
    ArrayView<iter> new_view{std::move(tmp_view)};

    REQUIRE(new_view.size() == view.size());
    REQUIRE(new_view.front() == view.front());

    constexpr test_t new_elem{static_cast<test_t>(300)};
    new_view.front() = new_elem;
    REQUIRE(view.front() == new_elem);
    REQUIRE(data.front() == new_elem);
  }

  SECTION("Move assignment operator")
  {
    auto tmp_view = view;
    ArrayView<iter> new_view{};
    new_view = std::move(tmp_view);
    REQUIRE(new_view.size() == view.size());
    REQUIRE(new_view.front() == view.front());

    constexpr test_t new_elem{static_cast<test_t>(300)};
    new_view.front() = new_elem;
    REQUIRE(view.front() == new_elem);
    REQUIRE(data.front() == new_elem);
  }

#ifdef NDEBUG
  SECTION("Invalid range in constructor")
  {
    REQUIRE_THROWS_WITH(ArrayView<vec_iter>(data.end(), data.begin()), "Invalid view range.");
  }

  SECTION("Out of bounds access with at()")
  {
    ArrayView<vec_iter> _view(data.begin(), data.end());

    using difference_type = vec_iter::difference_type;
    difference_type const too_high{static_cast<difference_type>(data.size() + 10)};
    difference_type const too_low{static_cast<difference_type>(-10)};

    REQUIRE_THROWS_WITH(_view.at(too_high), "Out of bounds.");
    REQUIRE_THROWS_WITH(_view.at(too_low), "Negative indices are not supported.");
  }
#endif
}
