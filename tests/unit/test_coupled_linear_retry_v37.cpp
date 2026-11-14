#include "cfdx/physics/coupled_linear_retry_policy_v37.h"
#include <iostream>

int main() {
    using cfdx::physics::CoupledLinearRetryV37;
    struct Sample { bool affects_linear; bool failure; bool predictor; bool fail_fast; };
    const Sample samples[] = {
        {false,true,false,true},
        {true,true,false,false},
        {false,false,false,false},
        {false,true,true,false},
        {false,false,true,false},
        {true,true,true,false},
    };
    for (const auto& sample : samples) {
        const CoupledLinearRetryV37 policy{sample.affects_linear,
            sample.failure,sample.predictor};
        if (policy.fail_fast() != sample.fail_fast ||
            policy.retry_can_help() != (sample.failure && !sample.fail_fast)) {
            std::cerr << "V37 fail-closed coupled linear retry mismatch\n";
            return 1;
        }
    }
    return 0;
}
