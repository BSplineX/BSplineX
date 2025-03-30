#ifndef BSPLINEX_DEFINES_HPP
#define BSPLINEX_DEFINES_HPP

#include <cassert>
#include <limits>

template <typename T>
constexpr void debugassert(T exp, std::string const &msg)
{
  assert(((void)msg, exp)); // Preserve original behavior
}

template <typename T>
constexpr void releaseassert(T exp, std::string const &msg)
{
  if (!exp)
  {
    throw std::runtime_error(msg);
  }
}

#ifdef BSPLINEX_DEBUG_LOG_CALL
#include <cstdio>
#define DEBUG_LOG_CALL() std::puts(__PRETTY_FUNCTION__);
#else
#define DEBUG_LOG_CALL() ;
#endif

namespace bsplinex::constants
{

constexpr size_t DENSE_MAX_COL = 512;

template <typename T>
constexpr auto ZERO = static_cast<T>(0);

template <typename T>
constexpr auto ONE = static_cast<T>(1);

template <typename T>
constexpr auto TOL_MULTIPLIER = static_cast<T>(100000);

template <typename T>
constexpr auto RTOL = std::numeric_limits<T>::epsilon() * TOL_MULTIPLIER<T>;

template <typename T>
constexpr auto ATOL = std::numeric_limits<T>::epsilon() * TOL_MULTIPLIER<T>;

} // namespace bsplinex::constants

#endif
