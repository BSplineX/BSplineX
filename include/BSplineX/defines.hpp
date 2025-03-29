#ifndef BSPLINEX_DEFINES_HPP
#define BSPLINEX_DEFINES_HPP

#include <cassert>

#define debugassert(exp, msg) assert(((void)msg, exp))
#define releaseassert(exp, msg)                                                                    \
  if (!(exp))                                                                                      \
    throw std::runtime_error(msg);

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

} // namespace bsplinex::constants

#endif
