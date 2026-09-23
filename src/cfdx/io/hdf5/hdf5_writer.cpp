// M0.9-T01 — HDF5 writer

#include "hdf5_writer.h"
#include "schema.h"

#include <H5public.h>
#include <H5Epublic.h>
#include <H5Dpublic.h>
#include <H5Fpublic.h>
#include <H5Gpublic.h>
#include <H5Tpublic.h>
#include <H5Spublic.h>
#include <H5Apublic.h>
#include <H5Ppublic.h>
#include <H5FDpublic.h>
#ifdef CFDX_ENABLE_PARALLEL_HDF5
#include <H5FDmpio.h>
#endif
#ifdef CFDX_ENABLE_PARALLEL_HDF5
#include <mpi.h>
#endif
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <string>
#include <cstdint>
#include <iomanip>
#include <sstream>

namespace cfdx {
namespace io {

static herr_t write_attr_str(hid_t loc_id, const char* name, const std::string& value);

static hid_t create_file_access_plist()
{
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    if (plist < 0) return -1;
#ifdef CFDX_ENABLE_PARALLEL_HDF5
    if (H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL) < 0) {
        H5Pclose(plist);
        return -1;
    }
#endif
    return plist;
}

static hid_t open_file_access_plist()
{
    return create_file_access_plist();
}

static std::uint64_t fnv1a_update(
    std::uint64_t hash, const void* data, std::size_t size)
{
    const auto* bytes = static_cast<const std::uint8_t*>(data);
    for (std::size_t i = 0; i < size; ++i) {
        hash ^= bytes[i];
        hash *= 1099511628211ULL;
    }
    return hash;
}

template<class T>
static std::uint64_t fnv1a_update_vector(
    std::uint64_t hash, const T* data, std::size_t count)
{
    // Empty CFDX containers may expose a null data pointer even when their
    // logical CSR offset array has the single zero entry. Do not dereference
    // such a pointer while constructing schema hashes.
    if (count == 0 || data == nullptr)
        return hash;
    return fnv1a_update(hash, data, count * sizeof(T));
}

static std::string hash_hex(std::uint64_t hash)
{
    std::ostringstream os;
    os << std::hex << std::setw(16) << std::setfill('0') << hash;
    return os.str();
}

static void write_schema_attributes(hid_t file, const cfdx::core::Mesh& mesh)
{
    write_attr_str(file, "format_version", std::to_string(CFDX_HDF5_FORMAT_VERSION));
    write_attr_str(file, "schema_version", std::to_string(CFDX_HDF5_SCHEMA_VERSION));
    write_attr_str(file, "cfdx_version", CFDX_VERSION);

    std::uint64_t topology = 1469598103934665603ULL;
    topology = fnv1a_update_vector(topology, mesh.faces().vertices_data(), mesh.faces().n_vertices());
    topology = fnv1a_update_vector(topology, mesh.faces().offsets_data(), mesh.faces().n_faces() + 1);
    topology = fnv1a_update_vector(topology, mesh.ownership().owner_data(), mesh.ownership().size());
    topology = fnv1a_update_vector(topology, mesh.ownership().neighbour_data(), mesh.ownership().size());
    topology = fnv1a_update_vector(topology, mesh.cells().faces_data(), mesh.cells().n_face_refs());
    topology = fnv1a_update_vector(topology, mesh.cells().offsets_data(), mesh.cells().n_cells() + 1);
    write_attr_str(file, "topology_hash", hash_hex(topology));

    std::uint64_t geometry = topology;
    geometry = fnv1a_update_vector(geometry, mesh.points().x_data(), mesh.n_points());
    geometry = fnv1a_update_vector(geometry, mesh.points().y_data(), mesh.n_points());
    geometry = fnv1a_update_vector(geometry, mesh.points().z_data(), mesh.n_points());
    write_attr_str(file, "mesh_hash", hash_hex(geometry));
}

static herr_t write_dataset(hid_t loc_id, const char* name,
                            const void* data, const hsize_t* dims, int rank) {
    hsize_t total = 1;
    for (int i = 0; i < rank; ++i) total *= dims[i];
    if (total == 0) {
        // For empty datasets, create with proper type but don't write data
        hid_t space = H5Screate_simple(rank, dims, nullptr);
        if (space < 0) return -1;
        hid_t ds = H5Dcreate2(loc_id, name, H5T_NATIVE_DOUBLE, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (ds < 0) { H5Sclose(space); return -1; }
        H5Dclose(ds);
        H5Sclose(space);
        return 0;
    }

    hid_t space = H5Screate_simple(rank, dims, nullptr);
    if (space < 0) return -1;
    hid_t ds = H5Dcreate2(loc_id, name, H5T_NATIVE_DOUBLE, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) { H5Sclose(space); return -1; }
    herr_t status = H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                             H5P_DEFAULT, data);
    H5Dclose(ds);
    H5Sclose(space);
    return status;
}

static herr_t write_dataset_u64(hid_t file_id, const char* name,
                                const std::uint64_t* data, hsize_t n) {
    if (n == 0) {
        hid_t space = H5Screate_simple(1, &n, nullptr);
        if (space < 0) return -1;
        hid_t type = H5T_NATIVE_UINT64;
        hid_t ds = H5Dcreate2(file_id, name, type, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (ds < 0) { H5Sclose(space); return -1; }
        H5Dclose(ds);
        H5Sclose(space);
        return 0;
    }

    hid_t space = H5Screate_simple(1, &n, nullptr);
    if (space < 0) return -1;
    hid_t type = H5T_NATIVE_UINT64;
    hid_t ds = H5Dcreate2(file_id, name, type, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) { H5Sclose(space); return -1; }
    herr_t status = H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(ds);
    H5Sclose(space);
    return status;
}

static herr_t write_dataset_i64(hid_t file_id, const char* name,
                                const std::int64_t* data, hsize_t n) {
    if (n == 0) {
        hid_t space = H5Screate_simple(1, &n, nullptr);
        if (space < 0) return -1;
        hid_t type = H5T_NATIVE_INT64;
        hid_t ds = H5Dcreate2(file_id, name, type, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (ds < 0) { H5Sclose(space); return -1; }
        H5Dclose(ds);
        H5Sclose(space);
        return 0;
    }

    hid_t space = H5Screate_simple(1, &n, nullptr);
    if (space < 0) return -1;
    hid_t type = H5T_NATIVE_INT64;
    hid_t ds = H5Dcreate2(file_id, name, type, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) { H5Sclose(space); return -1; }
    herr_t status = H5Dwrite(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(ds);
    H5Sclose(space);
    return status;
}

static herr_t write_attr_str(hid_t loc_id, const char* name, const std::string& value) {
    hid_t space = H5Screate(H5S_SCALAR);
    if (space < 0) return -1;
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, value.size() + 1);
    hid_t attr = H5Acreate2(loc_id, name, atype, space,
                            H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) { H5Sclose(space); H5Tclose(atype); return -1; }
    herr_t status = H5Awrite(attr, atype, value.c_str());
    H5Aclose(attr);
    H5Sclose(space);
    H5Tclose(atype);
    return status;
}

static bool ensure_group(hid_t file_id, const std::string& path) {
    std::vector<std::string> parts;
    size_t start = 0;
    if (!path.empty() && path[0] == '/') start = 1;
    
    while (start < path.size()) {
        size_t next = path.find('/', start);
        if (next == std::string::npos) next = path.size();
        
        if (next > start) {
            parts.emplace_back(path.substr(start, next - start));
        }
        start = next + 1;
    }
    
    std::string current_path;
    for (const auto& part : parts) {
        if (!current_path.empty()) current_path += "/";
        current_path += part;
        
        hid_t grp = H5Gopen2(file_id, current_path.c_str(), H5P_DEFAULT);
        if (grp < 0) {
            grp = H5Gcreate2(file_id, current_path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (grp < 0) return false;
        }
        if (grp >= 0) H5Gclose(grp);
    }
    return true;
}

bool write_mesh_hdf5(const std::string& filename, const cfdx::core::Mesh& mesh) {
    hid_t fapl = create_file_access_plist();
    if (fapl < 0) return false;
    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    H5Pclose(fapl);
    if (file < 0) return false;

    write_schema_attributes(file, mesh);

    write_attr_str(file, "n_points", std::to_string(mesh.n_points()));
    write_attr_str(file, "n_faces", std::to_string(mesh.n_faces()));
    write_attr_str(file, "n_cells", std::to_string(mesh.n_cells()));

    // Points (n_points x 3).
    {
        const hsize_t dims[2] = {mesh.n_points(), 3};
        std::vector<double> pts(mesh.n_points() * 3);
        const cfdx::core::PointCloud& pc = mesh.points();
        for (std::size_t i = 0; i < mesh.n_points(); ++i) {
            pts[i * 3 + 0] = pc.x(i);
            pts[i * 3 + 1] = pc.y(i);
            pts[i * 3 + 2] = pc.z(i);
        }
        write_dataset(file, "points", pts.data(), dims, 2);
    }

    // Faces (CSR) : vertices + offsets.
    {
        const cfdx::core::FaceConnectivity& fc = mesh.faces();
        write_dataset_u64(file, "face_vertices", fc.vertices_data(), fc.n_vertices());
        write_dataset_u64(file, "face_offsets", fc.offsets_data(), fc.n_faces() + 1);
    }

    // Owner / neighbour.
    {
        const cfdx::core::FaceOwnership& own = mesh.ownership();
        write_dataset_u64(file, "owner", own.owner_data(), own.size());
        write_dataset_i64(file, "neighbour", own.neighbour_data(), own.size());
    }

    // Cell faces (CSR).
    {
        const cfdx::core::CellConnectivity& cc = mesh.cells();
        write_dataset_u64(file, "cell_faces", cc.faces_data(), cc.n_face_refs());
        write_dataset_u64(file, "cell_offsets", cc.offsets_data(), cc.n_cells() + 1);
    }

    // Boundary patches.
    {
        const cfdx::core::BoundaryPatches& bp = mesh.boundary();
        std::size_t n_patches = bp.n_patches();
        std::vector<std::uint64_t> patch_starts(n_patches);
        std::vector<std::uint64_t> patch_counts(n_patches);
        std::vector<std::string> patch_names(n_patches);
        std::vector<std::uint32_t> patch_types(n_patches);
        
        for (std::size_t i = 0; i < n_patches; ++i) {
            const auto& patch = bp.patch(i);
            patch_starts[i] = patch.face_ids.empty() ? 0 : static_cast<std::uint64_t>(patch.face_ids[0]);
            patch_counts[i] = patch.face_ids.size();
            patch_names[i] = patch.name;
            patch_types[i] = static_cast<std::uint32_t>(patch.type);
        }
        
        // Boundary patch metadata is optional. For a mesh with no patches,
        // omit the attribute and datasets entirely; this keeps an empty mesh
        // representable without manufacturing an otherwise meaningless
        // boundary-patch schema.
        if (n_patches > 0) {
            std::string patches_str;
            for (std::size_t i = 0; i < n_patches; ++i) {
                if (i > 0) patches_str += ";";
                patches_str += patch_names[i] + ":" + std::to_string(patch_starts[i]) + ":" +
                               std::to_string(patch_counts[i]) + ":" + std::to_string(patch_types[i]);
            }
            write_attr_str(file, "boundary_patches", patches_str);

            std::vector<std::uint64_t> all_face_ids;
            for (std::size_t i = 0; i < n_patches; ++i) {
                const auto& patch = bp.patch(i);
                all_face_ids.insert(all_face_ids.end(), patch.face_ids.begin(), patch.face_ids.end());
            }
            std::vector<std::uint64_t> patch_face_offsets(n_patches + 1);
            patch_face_offsets[0] = 0;
            for (std::size_t i = 0; i < n_patches; ++i) {
                patch_face_offsets[i + 1] = patch_face_offsets[i] + patch_counts[i];
            }
            if (write_dataset_u64(file, "patch_face_ids", all_face_ids.data(),
                                  all_face_ids.size()) < 0 ||
                write_dataset_u64(file, "patch_face_offsets", patch_face_offsets.data(),
                                  patch_face_offsets.size()) < 0) {
                H5Fclose(file);
                return false;
            }
        }
    }

    H5Fclose(file);
    return true;
}

bool write_field_hdf5(const std::string& filename,
                       const cfdx::core::Field<double, cfdx::core::Location::CELL>& field) {
    hid_t fapl = open_file_access_plist();
    if (fapl < 0) return false;
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);
    H5Pclose(fapl);
    if (file < 0) return false;

    // Fixed group path: "fields"
    if (!ensure_group(file, "fields")) {
        H5Fclose(file);
        return false;
    }

    hid_t grp = H5Gopen2(file, "fields", H5P_DEFAULT);
    if (grp < 0) { H5Fclose(file); return false; }

    write_attr_str(grp, "name", field.name());
    write_attr_str(grp, "unit", field.metadata().unit);
    write_attr_str(grp, "dimension", std::to_string(field.dimension()));

    // Field uses SoA layout: need to interleave components for flat storage
    const std::size_t n = field.size();
    const std::size_t dim = field.dimension();
    std::vector<double> flat_values(n * dim);
    for (std::size_t c = 0; c < dim; ++c) {
        const double* comp_data = field.component_data(c);
        for (std::size_t i = 0; i < n; ++i) {
            flat_values[i * dim + c] = comp_data[i];
        }
    }

    const hsize_t dims[1] = {flat_values.size()};
    write_dataset(grp, "values", flat_values.data(), dims, 1);

    H5Gclose(grp);
    H5Fclose(file);
    return true;
}

}  // namespace io
}  // namespace cfdx
