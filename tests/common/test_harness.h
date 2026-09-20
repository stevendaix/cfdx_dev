// Test harness for CFDX unit tests
#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <functional>

namespace cfdx {
namespace testing {

struct TestCase {
    std::string name;
    std::function<void()> func;
};

std::vector<TestCase>& get_test_cases() {
    static std::vector<TestCase> cases;
    return cases;
}

void run_case(const std::string& name, std::function<void()> func) {
    get_test_cases().push_back({name, func});
}

#define EXPECT_TRUE(cond) \
    do { \
        if (!(cond)) { \
            std::cerr << "  FAILED: " << #cond << " in " << __FUNCTION__ << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            throw std::runtime_error("Assertion failed: " #cond); \
        } \
    } while (0)

#define EXPECT_FALSE(cond) \
    do { \
        if (cond) { \
            std::cerr << "  FAILED: " << #cond << " in " << __FUNCTION__ << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            throw std::runtime_error("Assertion failed: " #cond); \
        } \
    } while (0)

#define EXPECT_NEAR(a, b, eps) \
    do { \
        if (std::abs((a) - (b)) > (eps)) { \
            std::cerr << "  FAILED: " << #a " ≈ " #b " (diff=" << std::abs((a) - (b)) << " > " << eps << ")" << " in " << __FUNCTION__ << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            throw std::runtime_error("Assertion failed: " #a " ≈ " #b); \
        } \
    } while (0)

#define EXPECT_THROW(expr, exception_type) \
    do { \
        bool threw = false; \
        try { \
            expr; \
        } catch (const exception_type&) { \
            threw = true; \
        } catch (...) { \
            std::cerr << "  FAILED: wrong exception type thrown in " << __FUNCTION__ << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            throw std::runtime_error("Wrong exception type"); \
        } \
        if (!threw) { \
            std::cerr << "  FAILED: expected exception " #exception_type " not thrown in " << __FUNCTION__ << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            throw std::runtime_error("Expected exception not thrown"); \
        } \
    } while (0)

int run_all() {
    auto& cases = get_test_cases();
    int passed = 0;
    int failed = 0;

    std::cout << "Running " << cases.size() << " test cases..." << std::endl;

    for (const auto& tc : cases) {
        try {
            tc.func();
            std::cout << "  PASS: " << tc.name << std::endl;
            ++passed;
        } catch (const std::exception& e) {
            std::cout << "  FAIL: " << tc.name << " - " << e.what() << std::endl;
            ++failed;
        } catch (...) {
            std::cout << "  FAIL: " << tc.name << " - unknown exception" << std::endl;
            ++failed;
        }
    }

    std::cout << "\nResults: " << passed << " passed, " << failed << " failed" << std::endl;
    return failed == 0 ? 0 : 1;
}

}  // namespace testing
}  // namespace cfdx