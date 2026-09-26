#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/io/hdf5/hdf5_writer.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::testing;

int main() {
    run_case("hdf5_api_fails_explicitly_when_support_is_unavailable", [] {
        Mesh mesh;
        Field<double, Location::CELL> field(0, "value");
        EXPECT_FALSE(read_mesh_hdf5("unused.h5", mesh));
        EXPECT_FALSE(read_field_hdf5("unused.h5", field));
        EXPECT_FALSE(write_mesh_hdf5("unused.h5", mesh));
        EXPECT_FALSE(write_field_hdf5("unused.h5", field));
    });
    return run_all();
}
