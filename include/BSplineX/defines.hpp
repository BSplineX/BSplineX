#ifndef DEFINES_HPP
#define DEFINES_HPP

#include <cassert>

#define assertm(exp, msg) assert(((void)msg, exp))

#ifdef BSPLINEX_DEBUG_LOG_CALL
#include <cstdio>
#define DEBUG_LOG_CALL() std::puts(__PRETTY_FUNCTION__);
#else
#define DEBUG_LOG_CALL() ;
#endif

#define DENSE_MAX_COL 512

#endif
