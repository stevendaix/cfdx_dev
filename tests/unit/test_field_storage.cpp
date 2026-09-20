// M0.4-T02/T03 — Tests for FieldMetadata and StorageHandle

#include "cfdx/core/field/field.h"
#include "cfdx/core/field/storage.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    // --- Field metadata ---
    run_case("field_metadata_default", []() {
        ScalarCellField f(3, "p", "Pa", 1);
        EXPECT_TRUE(f.name() == "p");
        EXPECT_TRUE(f.metadata().unit == "Pa");
        EXPECT_TRUE(f.metadata().dimension == 1);
        EXPECT_TRUE(f.loc() == Location::CELL);
    });

    run_case("field_metadata_modifiable", []() {
        ScalarCellField f(3);
        f.metadata().name = "T";
        f.metadata().unit = "K";
        EXPECT_TRUE(f.name() == "T");
        EXPECT_TRUE(f.metadata().unit == "K");
    });

    // --- Storage states ---
    run_case("storage_default_host_only", []() {
        StorageHandle h;
        EXPECT_TRUE(h.state() == StorageState::HOST_ONLY);
        EXPECT_TRUE(h.size() == 0);
    });

    run_case("storage_resize", []() {
        StorageHandle h;
        h.resize(10);
        EXPECT_TRUE(h.size() == 10);
        EXPECT_FALSE(h.empty());
    });

    run_case("storage_clear", []() {
        StorageHandle h;
        h.resize(5);
        h.clear();
        EXPECT_TRUE(h.size() == 0);
        EXPECT_TRUE(h.empty());
        EXPECT_TRUE(h.state() == StorageState::HOST_ONLY);
    });

    run_case("storage_mark_host_dirty", []() {
        StorageHandle h;
        h.set_state(StorageState::SYNCHRONIZED);
        h.mark_host_dirty();
        EXPECT_TRUE(h.state() == StorageState::HOST_DIRTY);
    });

    run_case("storage_mark_device_dirty", []() {
        StorageHandle h;
        h.set_state(StorageState::SYNCHRONIZED);
        h.mark_device_dirty();
        EXPECT_TRUE(h.state() == StorageState::DEVICE_DIRTY);
    });

    run_case("storage_to_string", []() {
        EXPECT_TRUE(std::string(to_string(StorageState::HOST_ONLY)) == "host_only");
        EXPECT_TRUE(std::string(to_string(StorageState::DEVICE_ONLY)) == "device_only");
        EXPECT_TRUE(std::string(to_string(StorageState::SYNCHRONIZED)) == "synchronized");
        EXPECT_TRUE(std::string(to_string(StorageState::HOST_DIRTY)) == "host_dirty");
        EXPECT_TRUE(std::string(to_string(StorageState::DEVICE_DIRTY)) == "device_dirty");
    });

    run_case("location_to_string", []() {
        EXPECT_TRUE(std::string(to_string(Location::CELL)) == "cell");
        EXPECT_TRUE(std::string(to_string(Location::FACE)) == "face");
        EXPECT_TRUE(std::string(to_string(Location::POINT)) == "point");
        EXPECT_TRUE(std::string(to_string(Location::BOUNDARY)) == "boundary");
    });

    run_case("host_storage_data_access", []() {
        StorageHandle h;
        h.resize(3);
        h.host().data_ptr()[0] = 1.0;
        h.host().data_ptr()[1] = 2.0;
        h.host().data_ptr()[2] = 3.0;
        EXPECT_TRUE(h.host().data_ptr()[0] == 1.0);
        EXPECT_TRUE(h.host().data_ptr()[2] == 3.0);
    });

    return run_all();
}