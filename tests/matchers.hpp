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
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// BSplineX includes
#include "BSplineX/defines.hpp"

using namespace bsplinex::constants;

template <typename T>
class WithinAbsRelMatcher : public Catch::Matchers::MatcherBase<T>
{
public:
  WithinAbsRelMatcher(T const target, T const rtol, T const atol)
      : target(target), rtol(rtol), atol(atol),
        rel_matcher(Catch::Matchers::WithinRel(target, rtol)),
        abs_matcher(Catch::Matchers::WithinAbs(target, atol)), matcher(rel_matcher or abs_matcher)
  {
  }

  bool match(T const &actual) const override { return matcher.match(actual); }

  [[nodiscard]] std::string describe() const override { return matcher.describe(); }

private:
  T target;
  T rtol;
  T atol;
  Catch::Matchers::WithinRelMatcher rel_matcher;
  Catch::Matchers::WithinAbsMatcher abs_matcher;
  Catch::Matchers::Detail::MatchAnyOf<T> matcher;
};

template <typename T>
WithinAbsRelMatcher<T> WithinAbsRel(T expected, T rtol = RTOL<T>, T atol = ATOL<T>)
{
  return WithinAbsRelMatcher<T>(expected, rtol, atol);
}

template <typename T>
struct WithinAbsRelVectorMatcher final : Catch::Matchers::MatcherBase<std::vector<T>>
{
  static_assert(
      std::is_floating_point_v<T>,
      "WithinRelVectorMatcher can only be used with floating-point types"
  );

  WithinAbsRelVectorMatcher(std::vector<T> const &target, T const rtol, T const atol)
      : target(target), rtol(rtol), atol(atol)
  {
  }

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
      if (not Catch::Matchers::WithinRel(this->target[i], this->rtol).match(actual[i]) and
          not Catch::Matchers::WithinAbs(this->target[i], this->atol).match(actual[i]))
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
    ss << "is within " << this->rtol << " relative tolerance and " << this->atol
       << " absolute tolerance of expected vector";

    if (this->failures.empty() and this->mismatch_description.empty())
    {
      return "";
    }
    ss << "\n";

    if (not this->mismatch_description.empty())
    {
      ss << this->mismatch_description;
      return ss.str();
    }

    ss << "Vector comparison failed at " << this->failures.size() << " position(s):\n";

    for (auto const &failure : this->failures)
    {
      T relDiff = ZERO<T>;
      T absDiff = std::abs(failure.actual - failure.expected);
      if (failure.expected != ZERO<T>)
      {
        relDiff = std::abs((failure.actual - failure.expected) / failure.expected);
      }
      else
      {
        relDiff = std::abs(failure.actual) > std::numeric_limits<T>::epsilon()
                      ? std::numeric_limits<T>::infinity()
                      : ZERO<T>;
      }

      ss << "  [" << failure.index << "] Expected: " << failure.expected
         << ", Actual: " << failure.actual << ", Relative diff: " << relDiff << " (exceeds "
         << this->rtol << "), Absolute diff: " << absDiff << " (exceeds " << this->atol << ")\n";
    }

    return ss.str();
  }

private:
  std::vector<T> target;
  T rtol;
  T atol;

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
WithinAbsRelVectorMatcher<T>
VectorsWithinAbsRel(std::vector<T> const &expected, T rtol = RTOL<T>, T atol = ATOL<T>)
{
  return WithinAbsRelVectorMatcher<T>(expected, rtol, atol);
}

#endif // BSPLINEX_MATCHERS_HPP
