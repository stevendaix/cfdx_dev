#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/io/hdf5/hdf5_writer.h"
#include "cfdx/core/mesh/mesh.h"
#include "common/test_harness.h"

#include <hdf5.h>
#include <cstdio>
#include <string>

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::testing;

static Mesh cube() {
    Mesh m;
    m.points().resize(8);
    m.points().set(0,0,0,0); m.points().set(1,1,0,0);
    m.points().set(2,1,1,0); m.points().set(3,0,1,0);
    m.points().set(4,0,0,1); m.points().set(5,1,0,1);
    m.points().set(6,1,1,1); m.points().set(7,0,1,1);
    m.faces().push_face({0,3,2,1}); m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4}); m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3}); m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for(std::size_t f=0; f<6; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    return m;
}

static void replace_attr(hid_t file, const char* name, const char* value) {
    EXPECT_TRUE(H5Aexists(file, name) > 0);
    hid_t attr = H5Aopen(file, name, H5P_DEFAULT);
    EXPECT_TRUE(attr >= 0);
    if (attr < 0) return;
    hid_t type = H5Aget_type(attr);
    EXPECT_TRUE(type >= 0);
    if (type >= 0) EXPECT_TRUE(H5Awrite(attr, type, value) >= 0);
    if (type >= 0) H5Tclose(type);
    H5Aclose(attr);
}

static void write_and_close(const std::string& f) {
    std::remove(f.c_str());
    EXPECT_TRUE(write_mesh_hdf5(f, cube()));
}

int main() {
    run_case("reader_rejects_unsupported_schema_version", [] {
        const std::string f="/tmp/cfdx_schema_version.h5";
        write_and_close(f);
        hid_t file=H5Fopen(f.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
        EXPECT_TRUE(file>=0);
        if(file>=0) { replace_attr(file,"schema_version","999"); H5Fclose(file); }
        Mesh loaded;
        EXPECT_FALSE(read_mesh_hdf5(f,loaded));
        std::remove(f.c_str());
    });

    run_case("reader_rejects_format_version_change", [] {
        const std::string f="/tmp/cfdx_format_version.h5";
        write_and_close(f);
        hid_t file=H5Fopen(f.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
        EXPECT_TRUE(file>=0);
        if(file>=0) { replace_attr(file,"format_version","999"); H5Fclose(file); }
        Mesh loaded;
        EXPECT_FALSE(read_mesh_hdf5(f,loaded));
        std::remove(f.c_str());
    });

    run_case("reader_rejects_topology_hash_tampering", [] {
        const std::string f="/tmp/cfdx_topology_hash.h5";
        write_and_close(f);
        hid_t file=H5Fopen(f.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
        EXPECT_TRUE(file>=0);
        if(file>=0) { replace_attr(file,"topology_hash","0000000000000000"); H5Fclose(file); }
        Mesh loaded;
        EXPECT_FALSE(read_mesh_hdf5(f,loaded));
        std::remove(f.c_str());
    });

    run_case("reader_rejects_mesh_hash_tampering", [] {
        const std::string f="/tmp/cfdx_mesh_hash.h5";
        write_and_close(f);
        hid_t file=H5Fopen(f.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
        EXPECT_TRUE(file>=0);
        if(file>=0) { replace_attr(file,"mesh_hash","0000000000000000"); H5Fclose(file); }
        Mesh loaded;
        EXPECT_FALSE(read_mesh_hdf5(f,loaded));
        std::remove(f.c_str());
    });

    run_case("reader_rejects_missing_integrity_metadata", [] {
        const std::string f="/tmp/cfdx_missing_integrity.h5";
        write_and_close(f);
        hid_t file=H5Fopen(f.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
        EXPECT_TRUE(file>=0);
        if(file>=0) {
            EXPECT_TRUE(H5Adelete(file,"mesh_hash")>=0);
            H5Fclose(file);
        }
        Mesh loaded;
        EXPECT_FALSE(read_mesh_hdf5(f,loaded));
        std::remove(f.c_str());
    });

    return run_all();
}
