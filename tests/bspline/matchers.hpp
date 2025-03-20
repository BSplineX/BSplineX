#ifndef BSPLINEX_MATCHERS_HPP
#define BSPLINEX_MATCHERS_HPP
// Standard includes
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

// Third-party includes
#include <catch2/matchers/catch_matchers_floating_point.hpp>

template <typename T>
struct WithinRelVectorMatcher final : Catch::Matchers::MatcherBase<std::vector<T>>
{
  static_assert(
      std::is_floating_point_v<T>,
      "WithinRelVectorMatcher can only be used with floating-point types"
  );

  WithinRelVectorMatcher(std::vector<T> const &target, T eps) : target(target), eps(eps) {}

  bool match(std::vector<T> const &actual) const override
  {
    if (actual.size() != this->target.size())
    {
      this->mismatch_description = "Vector sizes don't match. "
                                   "Expected: " +
                                   std::to_string(this->target.size()) +
                                   ", Actual: " + std::to_string(actual.size());
      return false;
    }

    bool allMatch = true;
    this->failures.clear();

    for (size_t i = 0; i < actual.size(); ++i)
    {
      if (!Catch::Matchers::WithinRel(this->target[i], this->eps).match(actual[i]))
      {
        allMatch = false;
        this->failures.push_back({i, this->target[i], actual[i]});
      }
    }

    return allMatch;
  }

  std::string describe() const override
  {
    std::ostringstream ss;
    ss << "is within " << this->eps << " relative tolerance of expected vector";

    if (this->failures.empty() && this->mismatch_description.empty())
    {
      return "";
    }
    ss << "\n";

    if (!this->mismatch_description.empty())
    {
      ss << this->mismatch_description;
      return ss.str();
    }

    ss << "Vector comparison failed at " << this->failures.size() << " position(s):\n";

    for (auto const &failure : this->failures)
    {
      T relDiff = static_cast<T>(0);
      if (failure.expected != static_cast<T>(0))
      {
        relDiff = std::abs((failure.actual - failure.expected) / failure.expected);
      }
      else
      {
        relDiff = std::abs(failure.actual) > std::numeric_limits<T>::epsilon()
                      ? std::numeric_limits<T>::infinity()
                      : static_cast<T>(0);
      }

      ss << "  [" << failure.index << "] Expected: " << failure.expected
         << ", Actual: " << failure.actual << ", Relative diff: " << relDiff << " (exceeds " << eps
         << ")\n";
    }

    return ss.str();
  }

private:
  std::vector<T> target;
  T eps;

  struct Failure
  {
    size_t index;
    T expected;
    T actual;
  };

  mutable std::vector<Failure> failures;
  mutable std::string mismatch_description;
};

template <typename T>
WithinRelVectorMatcher<T>
VectorsWithinRel(std::vector<T> const &expected, T margin = std::numeric_limits<T>::epsilon() * 100)
{
  return WithinRelVectorMatcher<T>(expected, margin);
}

#endif // BSPLINEX_MATCHERS_HPP
