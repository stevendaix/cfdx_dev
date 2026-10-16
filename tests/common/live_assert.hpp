#pragma once
// Live assertion for Release builds.
//
// <cassert> compiles assert() to a no-op when NDEBUG is defined, which is
// the case in every Release/CI build of this repository
// (CMAKE_CXX_FLAGS_RELEASE and the GNU/Clang branch of the top-level
// CMakeLists both add -DNDEBUG). Several tests relied on assert() for
// their acceptance checks, so those checks were silently inert in
// Release. The 2026-10-16 audit (VMFL036 kernel test) showed the failure
// mode: a corrupted control value passed the Release test suite.
//
// This header redefines assert() into a check that always throws, so the
// tests are fail-closed in every build type. Include it AFTER the last
// #include of the translation unit so the macro is (re)defined after any
// inclusion of <cassert>.
#include <cassert>
#include <stdexcept>
#include <string>

#undef assert
#define assert(expr)                                                     \
    do {                                                                 \
        if (!(expr)) {                                                  \
            throw std::runtime_error(                                   \
                "live_assert failed: " #expr " (" __FILE__ ":" +         \
                std::to_string(__LINE__) + ")");                          \
        }                                                                \
    } while (0)
