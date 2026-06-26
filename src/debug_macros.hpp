#ifndef DEBUG_MACRO_HPP
#define DEBUG_MACRO_HPP

#include <atomic>
#include <chrono>
#ifndef __DEBUG
#define __DEBUG 1 
#endif

#if __DEBUG 
#define DEBUG(s) do { \
    std::cerr << "[" << __FILE__ << "][" << __FUNCTION__ << "][" << __LINE__ << "]: " \
              << s << std::endl; \
} while (0)
#else 
#define DEBUG(s)
#endif

#ifdef __DEBUG

#include <atomic>
#include <chrono>
#include <iostream>

inline std::atomic<std::chrono::high_resolution_clock::time_point> dbg_timer;

#define TIME_DEBUG(s) do { \
    auto _now = std::chrono::high_resolution_clock::now(); \
    double _ms = std::chrono::duration<double, std::milli>(_now - dbg_timer.load()).count(); \
    std::cerr << "[" << __FILE__ << "][" << __FUNCTION__ << "][" << __LINE__ \
              << "]: " << _ms << "ms - " << s << std::endl; \
    dbg_timer.store(_now); \
} while (0)

#define TIME_DEBUG_INIT() do { \
    dbg_timer.store(std::chrono::high_resolution_clock::now()); \
} while (0)

#else

#define TIME_DEBUG(s)      do {} while (0)
#define TIME_DEBUG_INIT()  do {} while (0)

#endif

#endif // !DEBUG_MACRO_HPP
