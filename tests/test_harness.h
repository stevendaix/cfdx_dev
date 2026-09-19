#pragma once

// CFDX minimal test harness — NOT a test framework project.
// Only the macros needed for CFDX tests. Integrated with CMake ctest.

#include <cstdio>
#include <string>
#include <functional>
#include <stdexcept>
#include <cmath>

namespace cfdx {
namespace testing {

struct TestFailure : std::exception {
    explicit TestFailure(std::string msg) : message(std::move(msg)) {}
    const char* what() const noexcept override { return message.c_str(); }
    std::string message;
};

inline int& run_count() { static int c = 0; return c; }
inline int& fail_count() { static int c = 0; return c; }

inline void run_case(const char* name, std::function<void()> fn) {
    ++run_count();
    std::fprintf(stdout, "[test] %s\n", name);
    try {
        fn();
    } catch (const TestFailure& e) {
        ++fail_count();
        std::fprintf(stderr, "  FAIL: %s\n", e.what());
    } catch (const std::exception& e) {
        ++fail_count();
        std::fprintf(stderr, "  FAIL (exception): %s\n", e.what());
    }
}

#define TEST(name, body) \
    cfdx::testing::run_case((name), []() -> void { body; })

inline void expect(bool cond, const char* expr, const char* file, int line) {
    if (!cond) {
        std::fprintf(stderr, "  assertion failed: %s at %s:%d\n", expr, file, line);
        throw TestFailure(std::string("assertion failed: ") + expr + " at " + file + ":" + std::to_string(line));
    }
}

#define EXPECT_TRUE(c)  cfdx::testing::expect((c), #c, __FILE__, __LINE__)
#define EXPECT_FALSE(c) cfdx::testing::expect(!(c), "!(" #c ")", __FILE__, __LINE__)

#define EXPECT_EQ(a, b)                                                    \
    do {                                                                   \
        auto _la = (a);                                                    \
        auto _lb = (b);                                                    \
        if (!(_la == _lb)) {                                               \
            std::fprintf(stderr, "  EXPECT_EQ failed: %s != %s at %s:%d\n", \
                         #a, #b, __FILE__, __LINE__);                      \
            throw cfdx::testing::TestFailure(                               \
                std::string("EXPECT_EQ: ") + #a + " != " + #b +            \
                " at " + __FILE__ + ":" + std::to_string(__LINE__));       \
        }                                                                  \
    } while (0)

#define EXPECT_NE(a, b)                                                    \
    do {                                                                   \
        auto _la = (a);                                                    \
        auto _lb = (b);                                                    \
        if (_la == _lb) {                                                  \
            std::fprintf(stderr, "  EXPECT_NE failed: %s == %s at %s:%d\n", \
                         #a, #b, __FILE__, __LINE__);                      \
            throw cfdx::testing::TestFailure(                               \
                std::string("EXPECT_NE: ") + #a + " == " + #b +            \
                " at " + __FILE__ + ":" + std::to_string(__LINE__));       \
        }                                                                  \
    } while (0)

#define EXPECT_NEAR(a, b, tol)                                             \
    do {                                                                   \
        auto _la = (a);                                                    \
        auto _lb = (b);                                                    \
        auto _diff = std::abs(static_cast<long double>(_la) -              \
                              static_cast<long double>(_lb));              \
        if (!(_diff <= static_cast<long double>(tol))) {                   \
            std::fprintf(stderr, "  EXPECT_NEAR failed: |%s - %s| = %Lg > %s " \
                         "at %s:%d\n",                                     \
                         #a, #b, _diff, #tol, __FILE__, __LINE__);         \
            throw cfdx::testing::TestFailure(                               \
                std::string("EXPECT_NEAR: |") + #a + " - " + #b +           \
                std::string("| > ") + #tol +                                \
                " at " + __FILE__ + ":" + std::to_string(__LINE__));       \
        }                                                                  \
    } while (0)

#define EXPECT_THROW(stmt, exc)                                            \
    do {                                                                   \
        bool _caught = false;                                              \
        try {                                                              \
            stmt;                                                          \
        } catch (const exc&) {                                             \
            _caught = true;                                                \
        } catch (...) {                                                    \
        }                                                                  \
        if (!_caught) {                                                    \
            std::fprintf(stderr, "  EXPECT_THROW failed: " #stmt           \
                         " did not throw " #exc " at %s:%d\n",            \
                         __FILE__, __LINE__);                              \
            throw cfdx::testing::TestFailure(                               \
                std::string("EXPECT_THROW: ") + #stmt +                    \
                " did not throw " + #exc +                                 \
                " at " + __FILE__ + ":" + std::to_string(__LINE__));       \
        }                                                                  \
    } while (0)

inline int run_all() {
    std::fprintf(stdout, "Tests: %d run, %d failed\n", run_count(), fail_count());
    return fail_count() == 0 ? 0 : 1;
}

}  // namespace testing
}  // namespace cfdx