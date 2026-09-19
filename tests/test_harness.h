#pragma once

#include <cstdio>
#include <string>
#include <functional>
#include <stdexcept>

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

inline void expect(bool cond, const char* expr, const char* file, int line) {
    if (!cond) {
        std::fprintf(stderr, "  assertion failed: %s at %s:%d\n", expr, file, line);
        throw TestFailure(std::string("assertion failed: ") + expr + " at " + file + ":" + std::to_string(line));
    }
}

#define EXPECT(cond) cfdx::testing::expect((cond), #cond, __FILE__, __LINE__)

inline int run_all() {
    std::fprintf(stdout, "Tests: %d run, %d failed\n", run_count(), fail_count());
    return fail_count() == 0 ? 0 : 1;
}

}  // namespace testing
}  // namespace cfdx