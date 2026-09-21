#include "cfdx/runtime/gpu/cuda_backend.h"
#include "cfdx/runtime/gpu/cuda_field_mirror.h"
#include "cfdx/core/field/field.h"
#include "common/test_harness.h"

using namespace cfdx::testing;
using namespace cfdx::core;
using namespace cfdx::runtime::gpu;

int main() {
    run_case("cuda_field_roundtrip", [] {
        ScalarCellField f(8, "phi", "1", 1);
        for (std::size_t i=0; i<f.size(); ++i) f(i)=static_cast<double>(i)+0.5;
        CudaFieldMirror<Location::CELL> mirror;
        mirror.upload(f);
        f.fill(0.0);
        mirror.download(f);
        EXPECT_NEAR(f(0),0.5,1e-14);
        EXPECT_NEAR(f(7),7.5,1e-14);
    });

    run_case("cuda_nonblocking_stream", [] {
        CudaStream stream;
        EXPECT_TRUE(stream.get()!=nullptr);
        stream.synchronize();
    });

    return run_all();
}
