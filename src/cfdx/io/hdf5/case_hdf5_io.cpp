// M0.10-T04: CFDX case-level HDF5 I/O implementation
#include "case_hdf5_io.h"
#include "schema.h"
#include "mini_json.h"

#include <H5public.h>
#include <H5Epublic.h>
#include <H5Dpublic.h>
#include <H5Fpublic.h>
#include <H5Gpublic.h>
#include <H5Tpublic.h>
#include <H5Spublic.h>
#include <H5Apublic.h>
#include <H5Ppublic.h>
#include <hdf5.h>
#ifdef CFDX_ENABLE_PARALLEL_HDF5
#include <H5FDmpi.h>
#include <H5FDmpio.h>
#include <mpi.h>
#endif
#include <cstdio>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <string>
#include <cstdint>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <iostream>

namespace cfdx {
namespace io {

// ---------------------------------------------------------------------------
// File-access helper (parallel-aware, mirrors hdf5_reader/writer)
// ---------------------------------------------------------------------------
static hid_t case_create_file_access_plist() {
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    if (plist < 0) return -1;
#ifdef CFDX_ENABLE_PARALLEL_HDF5
    int mpi_initialized = 0;
    int mpi_finalized = 0;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized) MPI_Finalized(&mpi_finalized);
    if (mpi_initialized && !mpi_finalized) {
        if (H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL) < 0) {
            H5Pclose(plist);
            return -1;
        }
    }
#endif
    return plist;
}

// ---------------------------------------------------------------------------
// Low-level dataset readers (mirrors hdf5_reader.cpp)
// ---------------------------------------------------------------------------
static bool case_read_dataset_double(hid_t loc_id, const char* name,
                                std::vector<double>& out, int expected_rank = -1) {
    hid_t ds = H5Dopen2(loc_id, name, H5P_DEFAULT);
    if (ds < 0) return false;
    hid_t space = H5Dget_space(ds);
    hid_t type = H5Dget_type(ds);
    if (space < 0 || type < 0) {
        if (type >= 0) H5Tclose(type);
        if (space >= 0) H5Sclose(space);
        H5Dclose(ds);
        return false;
    }
    const int rank = H5Sget_simple_extent_ndims(space);
    if (rank < 0 || (expected_rank >= 0 && rank != expected_rank) ||
        H5Tget_class(type) != H5T_FLOAT || H5Tget_size(type) != sizeof(double)) {
        H5Tclose(type); H5Sclose(space); H5Dclose(ds);
        return false;
    }
    std::vector<hsize_t> dims(static_cast<std::size_t>(rank));
    H5Sget_simple_extent_dims(space, dims.data(), nullptr);
    hsize_t total = 1;
    for (hsize_t d : dims) total *= d;
    out.resize(static_cast<std::size_t>(total));
    herr_t status = total == 0 ? 0 :
        H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
    H5Tclose(type); H5Sclose(space); H5Dclose(ds);
    return status >= 0;
}

static bool case_read_dataset_u64(hid_t loc_id, const char* name,
                             std::vector<std::uint64_t>& out) {
    hid_t ds = H5Dopen2(loc_id, name, H5P_DEFAULT);
    if (ds < 0) return false;
    hid_t space = H5Dget_space(ds);
    hid_t type = H5Dget_type(ds);
    if (space < 0 || type < 0) {
        if (type >= 0) H5Tclose(type);
        if (space >= 0) H5Sclose(space);
        H5Dclose(ds);
        return false;
    }
    const int rank = H5Sget_simple_extent_ndims(space);
    if (rank != 1 || H5Tget_class(type) != H5T_INTEGER ||
        H5Tget_size(type) != sizeof(std::uint64_t) ||
        H5Tget_sign(type) != H5T_SGN_NONE) {
        H5Tclose(type); H5Sclose(space); H5Dclose(ds);
        return false;
    }
    hsize_t dim = 0;
    H5Sget_simple_extent_dims(space, &dim, nullptr);
    out.resize(static_cast<std::size_t>(dim));
    herr_t status = dim == 0 ? 0 :
        H5Dread(ds, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
    H5Tclose(type); H5Sclose(space); H5Dclose(ds);
    return status >= 0;
}

static bool case_read_dataset_i64(hid_t loc_id, const char* name,
                             std::vector<std::int64_t>& out) {
    hid_t ds = H5Dopen2(loc_id, name, H5P_DEFAULT);
    if (ds < 0) return false;
    hid_t space = H5Dget_space(ds);
    hid_t type = H5Dget_type(ds);
    if (space < 0 || type < 0) {
        if (type >= 0) H5Tclose(type);
        if (space >= 0) H5Sclose(space);
        H5Dclose(ds);
        return false;
    }
    const int rank = H5Sget_simple_extent_ndims(space);
    if (rank != 1 || H5Tget_class(type) != H5T_INTEGER ||
        H5Tget_size(type) != sizeof(std::int64_t) ||
        H5Tget_sign(type) != H5T_SGN_2) {
        H5Tclose(type); H5Sclose(space); H5Dclose(ds);
        return false;
    }
    hsize_t dim = 0;
    H5Sget_simple_extent_dims(space, &dim, nullptr);
    out.resize(static_cast<std::size_t>(dim));
    herr_t status = dim == 0 ? 0 :
        H5Dread(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
    H5Tclose(type); H5Sclose(space); H5Dclose(ds);
    return status >= 0;
}

// ---------------------------------------------------------------------------
// Attribute reader
// ---------------------------------------------------------------------------
static bool case_read_attr_str(hid_t loc_id, const char* name, std::string& out) {
    if (H5Aexists(loc_id, name) <= 0) return false;
    hid_t attr = H5Aopen(loc_id, name, H5P_DEFAULT);
    if (attr < 0) return false;

    hid_t atype = H5Aget_type(attr);
    if (atype < 0) { H5Aclose(attr); return false; }

    hsize_t type_size = H5Tget_size(atype);
    if (type_size == 0) {
        type_size = H5Aget_storage_size(attr);
        if (type_size == 0) type_size = 1;
    }
    std::size_t buf_size = static_cast<std::size_t>(type_size) + 1;
    std::vector<char> buf(buf_size, '\0');

    herr_t status = H5Aread(attr, atype, buf.data());
    H5Tclose(atype);
    H5Aclose(attr);

    if (status < 0) return false;

    out = buf.data();
    out.erase(out.find_last_not_of('\0') + 1);
    if (out.empty()) out = std::string(buf.data(), 1);
    return true;
}

// ---------------------------------------------------------------------------
// FNV-1a hashing (mirrors hdf5_reader/writer.cpp)
// ---------------------------------------------------------------------------
static std::uint64_t case_fnv1a_update(std::uint64_t hash, const void* data, std::size_t size) {
    const auto* bytes = static_cast<const std::uint8_t*>(data);
    for (std::size_t i = 0; i < size; ++i) {
        hash ^= bytes[i];
        hash *= 1099511628211ULL;
    }
    return hash;
}

template<class T>
static std::uint64_t case_fnv1a_update_vector(std::uint64_t hash, const T* data, std::size_t count) {
    if (count == 0 || data == nullptr) return hash;
    return case_fnv1a_update(hash, data, count * sizeof(T));
}

static std::string case_hash_hex(std::uint64_t hash) {
    std::ostringstream os;
    os << std::hex << std::setw(16) << std::setfill('0') << hash;
    return os.str();
}

// ---------------------------------------------------------------------------
// Dataset writers (mirrors hdf5_writer.cpp)
// ---------------------------------------------------------------------------
static herr_t case_write_attr_str(hid_t loc_id, const char* name, const std::string& value) {
    hid_t space = H5Screate(H5S_SCALAR);
    if (space < 0) return -1;
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, value.size() + 1);
    hid_t attr = H5Acreate2(loc_id, name, atype, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) { H5Sclose(space); H5Tclose(atype); return -1; }
    herr_t status = H5Awrite(attr, atype, value.c_str());
    H5Aclose(attr);
    H5Sclose(space);
    H5Tclose(atype);
    return status;
}

static herr_t case_write_dataset_u64(hid_t file_id, const char* name,
                                const std::uint64_t* data, hsize_t n) {
    if (n == 0) {
        hsize_t zero = 0;
        hid_t space = H5Screate_simple(1, &zero, nullptr);
        if (space < 0) return -1;
        hid_t ds = H5Dcreate2(file_id, name, H5T_NATIVE_UINT64, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (ds < 0) { H5Sclose(space); return -1; }
        H5Dclose(ds); H5Sclose(space); return 0;
    }
    hid_t space = H5Screate_simple(1, &n, nullptr);
    if (space < 0) return -1;
    hid_t ds = H5Dcreate2(file_id, name, H5T_NATIVE_UINT64, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) { H5Sclose(space); return -1; }
    herr_t status = H5Dwrite(ds, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(ds); H5Sclose(space);
    return status;
}

static herr_t case_write_dataset_i64(hid_t file_id, const char* name,
                                const std::int64_t* data, hsize_t n) {
    if (n == 0) {
        hsize_t zero = 0;
        hid_t space = H5Screate_simple(1, &zero, nullptr);
        if (space < 0) return -1;
        hid_t ds = H5Dcreate2(file_id, name, H5T_NATIVE_INT64, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (ds < 0) { H5Sclose(space); return -1; }
        H5Dclose(ds); H5Sclose(space); return 0;
    }
    hid_t space = H5Screate_simple(1, &n, nullptr);
    if (space < 0) return -1;
    hid_t ds = H5Dcreate2(file_id, name, H5T_NATIVE_INT64, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) { H5Sclose(space); return -1; }
    herr_t status = H5Dwrite(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(ds); H5Sclose(space);
    return status;
}

static herr_t write_dataset_double_2d(hid_t file_id, const char* name,
                                      const double* data, hsize_t rows, hsize_t cols) {
    const hsize_t dims[2] = {rows, cols};
    hsize_t total = rows * cols;
    if (total == 0) {
        hid_t space = H5Screate_simple(2, dims, nullptr);
        if (space < 0) return -1;
        hid_t ds = H5Dcreate2(file_id, name, H5T_NATIVE_DOUBLE, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (ds < 0) { H5Sclose(space); return -1; }
        H5Dclose(ds); H5Sclose(space); return 0;
    }
    hid_t space = H5Screate_simple(2, dims, nullptr);
    if (space < 0) return -1;
    hid_t ds = H5Dcreate2(file_id, name, H5T_NATIVE_DOUBLE, space,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (ds < 0) { H5Sclose(space); return -1; }
    herr_t status = H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(ds); H5Sclose(space);
    return status;
}

// ---------------------------------------------------------------------------
// Patch metadata parsing
// ---------------------------------------------------------------------------
static std::vector<std::string> case_parse_patch_metadata(const std::string& patches_str) {
    std::vector<std::string> result;
    if (patches_str.empty()) return result;
    size_t start = 0;
    while (start < patches_str.size()) {
        size_t end = patches_str.find(';', start);
        if (end == std::string::npos) end = patches_str.size();
        result.emplace_back(patches_str.substr(start, end - start));
        start = end + 1;
    }
    return result;
}

// ===========================================================================
// CaseSetup JSON serialization helpers
// ===========================================================================

static void case_setup_to_json(const CaseSetup& setup, mini_json::value& obj) {
    obj.type = mini_json::value_type::object;

    obj.object["schema_version"] = mini_json::value(static_cast<double>(setup.schema_version));

    // SourceInfo
    mini_json::value source_obj(mini_json::value_type::object);
    source_obj.object["solver"] = setup.source.solver;
    source_obj.object["version"] = setup.source.version;
    source_obj.object["case_path"] = setup.source.case_path;
    source_obj.object["case_name"] = setup.source.case_name;
    source_obj.object["format"] = setup.source.format;
    obj.object["source"] = std::move(source_obj);

    // MeshMetadata
    mini_json::value mesh_obj(mini_json::value_type::object);
    mesh_obj.object["n_vertices"] = mini_json::value(static_cast<double>(setup.mesh_info.n_vertices));
    mesh_obj.object["n_faces"] = mini_json::value(static_cast<double>(setup.mesh_info.n_faces));
    mesh_obj.object["n_cells"] = mini_json::value(static_cast<double>(setup.mesh_info.n_cells));
    mesh_obj.object["n_boundary_faces"] = mini_json::value(static_cast<double>(setup.mesh_info.n_boundary_faces));
    mesh_obj.object["n_internal_faces"] = mini_json::value(static_cast<double>(setup.mesh_info.n_internal_faces));
    mesh_obj.object["n_patches"] = mini_json::value(static_cast<double>(setup.mesh_info.n_patches));
    mesh_obj.object["dimension"] = mini_json::value(static_cast<double>(setup.mesh_info.dimension));
    mesh_obj.object["has_polyhedral"] = mini_json::value(setup.mesh_info.has_polyhedral);
    mesh_obj.object["has_nonplanar_faces"] = mini_json::value(setup.mesh_info.has_nonplanar_faces);
    mesh_obj.object["mesh_type"] = setup.mesh_info.mesh_type;

    mini_json::array_t cell_types_arr;
    for (const auto& ct : setup.mesh_info.cell_types) {
        cell_types_arr.emplace_back(ct);
    }
    mesh_obj.object["cell_types"] = std::move(cell_types_arr);
    obj.object["mesh_info"] = std::move(mesh_obj);

    // Physics
    obj.object["physics_model"] = setup.physics_model;
    obj.object["turbulence_model"] = setup.turbulence_model;
    obj.object["energy_model"] = setup.energy_model;
    obj.object["multiphase_model"] = setup.multiphase_model;
    obj.object["radiation_model"] = setup.radiation_model;
    obj.object["transient"] = mini_json::value(setup.transient);
    obj.object["time_step"] = mini_json::value(setup.time_step);
    obj.object["end_time"] = mini_json::value(setup.end_time);
    obj.object["max_time_steps"] = mini_json::value(static_cast<double>(setup.max_time_steps));

    // Materials
    mini_json::array_t mats_arr;
    for (const auto& mat : setup.materials) {
        mini_json::value m_obj(mini_json::value_type::object);
        m_obj.object["name"] = mat.name;
        m_obj.object["density"] = mini_json::value(mat.density);
        m_obj.object["dynamic_viscosity"] = mini_json::value(mat.dynamic_viscosity);
        m_obj.object["thermal_conductivity"] = mini_json::value(mat.thermal_conductivity);
        m_obj.object["specific_heat"] = mini_json::value(mat.specific_heat);
        m_obj.object["molecular_weight"] = mini_json::value(mat.molecular_weight);
        m_obj.object["eos_model"] = mat.eos_model;
        m_obj.object["turbulence_model"] = mat.turbulence_model;

        mini_json::value extras_obj(mini_json::value_type::object);
        for (const auto& [k, v] : mat.extra_properties) {
            extras_obj.object[k] = v;
        }
        m_obj.object["extra_properties"] = std::move(extras_obj);
        mats_arr.push_back(std::move(m_obj));
    }
    obj.object["materials"] = std::move(mats_arr);

    // Boundary conditions
    mini_json::array_t bcs_arr;
    for (const auto& bc : setup.boundary_conditions) {
        mini_json::value bc_obj(mini_json::value_type::object);
        bc_obj.object["patch_name"] = bc.patch_name;

        auto bc_type_str = [](BCType t) -> const char* {
            switch (t) {
                case BCType::WALL:            return "wall";
                case BCType::INLET:           return "inlet";
                case BCType::OUTLET:          return "outlet";
                case BCType::PRESSURE_OUTLET: return "pressure_outlet";
                case BCType::SYMMETRY:        return "symmetry";
                case BCType::PERIODIC:        return "periodic";
                case BCType::INTERFACE:       return "interface";
                case BCType::EMPTY:           return "empty";
                case BCType::INTERNAL:        return "internal";
                default:                      return "unknown";
            }
        };
        auto bcvt_str = [](BCValueType t) -> const char* {
            switch (t) {
                case BCValueType::FIXED:            return "fixed";
                case BCValueType::ZERO_GRADIENT:    return "zero_gradient";
                case BCValueType::MIXED:            return "mixed";
                case BCValueType::OUTLET_PRESSURE:  return "outlet_pressure";
                case BCValueType::WALL_NO_SLIP:     return "wall_no_slip";
                case BCValueType::WALL_SLIP:        return "wall_slip";
                case BCValueType::WALL_THERMAL:     return "wall_thermal";
                default:                            return "unknown";
            }
        };

        bc_obj.object["type"] = bc_type_str(bc.type);
        bc_obj.object["value_type"] = bcvt_str(bc.value_type);
        bc_obj.object["velocity_magnitude"] = mini_json::value(bc.velocity_magnitude);

        mini_json::array_t vv_arr;
        for (double v : bc.velocity_vector) {
            vv_arr.emplace_back(mini_json::value(v));
        }
        bc_obj.object["velocity_vector"] = std::move(vv_arr);
        bc_obj.object["temperature"] = mini_json::value(bc.temperature);
        bc_obj.object["pressure"] = mini_json::value(bc.pressure);
        bc_obj.object["turbulence_intensity"] = mini_json::value(bc.turbulence_intensity);
        bc_obj.object["turbulence_length_scale"] = mini_json::value(bc.turbulence_length_scale);
        bc_obj.object["roughness_height"] = mini_json::value(bc.roughness_height);
        bc_obj.object["wall_temperature"] = mini_json::value(bc.wall_temperature);
        bc_obj.object["back_pressure"] = mini_json::value(bc.back_pressure);
        bc_obj.object["source_zone_name"] = bc.source_zone_name;
        bc_obj.object["source_type_name"] = bc.source_type_name;

        mini_json::value raw_obj(mini_json::value_type::object);
        for (const auto& [k, v] : bc.raw_params) {
            raw_obj.object[k] = v;
        }
        bc_obj.object["raw_params"] = std::move(raw_obj);
        bcs_arr.push_back(std::move(bc_obj));
    }
    obj.object["boundary_conditions"] = std::move(bcs_arr);

    // Initial conditions
    mini_json::value ic_obj(mini_json::value_type::object);
    ic_obj.object["velocity"] = mini_json::value(setup.initial_condition.velocity);
    ic_obj.object["pressure"] = mini_json::value(setup.initial_condition.pressure);
    ic_obj.object["temperature"] = mini_json::value(setup.initial_condition.temperature);

    mini_json::array_t ic_vv;
    for (double v : setup.initial_condition.velocity_vector) {
        ic_vv.emplace_back(mini_json::value(v));
    }
    ic_obj.object["velocity_vector"] = std::move(ic_vv);

    mini_json::value sf_obj(mini_json::value_type::object);
    for (const auto& [k, v] : setup.initial_condition.scalar_fields) {
        sf_obj.object[k] = mini_json::value(v);
    }
    ic_obj.object["scalar_fields"] = std::move(sf_obj);
    obj.object["initial_condition"] = std::move(ic_obj);

    // Numerics
    mini_json::value num_obj(mini_json::value_type::object);
    num_obj.object["momentum_scheme"] = setup.numerics.momentum_scheme;
    num_obj.object["pressure_scheme"] = setup.numerics.pressure_scheme;
    num_obj.object["momentum_interpolation"] = setup.numerics.momentum_interpolation;
    num_obj.object["transient_scheme"] = setup.numerics.transient_scheme;
    num_obj.object["gradient_operator"] = setup.numerics.gradient_operator;
    num_obj.object["under_relaxation_momentum"] = mini_json::value(setup.numerics.under_relaxation_momentum);
    num_obj.object["under_relaxation_pressure"] = mini_json::value(setup.numerics.under_relaxation_pressure);
    num_obj.object["coupled_solver"] = setup.numerics.coupled_solver;
    num_obj.object["preconditioner"] = setup.numerics.preconditioner;
    num_obj.object["residual_target"] = setup.numerics.residual_target;
    num_obj.object["max_iterations"] = mini_json::value(static_cast<double>(setup.numerics.max_iterations));

    mini_json::value raw_obj(mini_json::value_type::object);
    for (const auto& [k, v] : setup.numerics.raw_settings) {
        raw_obj.object[k] = v;
    }
    num_obj.object["raw_settings"] = std::move(raw_obj);
    // Canonical N1 numerical selections. This is deliberately separate from
    // legacy human-readable scheme fields: the registry keys are the machine
    // contract consumed by the numerical resolver.
    if (setup.has_explicit_numerics) {
        mini_json::value selection_obj(mini_json::value_type::object);
        mini_json::array_t entries_arr;
        for (const auto& entry : setup.numerical_config.entries) {
            mini_json::value entry_obj(mini_json::value_type::object);
            entry_obj.object["family"] = std::string(cfdx::core::to_string(entry.family));
            entry_obj.object["configuration_key"] = entry.configuration_key;
            entries_arr.push_back(std::move(entry_obj));
        }
        selection_obj.object["entries"] = std::move(entries_arr);
        mini_json::array_t required_arr;
        for (const auto family : setup.numerical_config.required_families)
            required_arr.emplace_back(std::string(cfdx::core::to_string(family)));
        selection_obj.object["required_families"] = std::move(required_arr);
        num_obj.object["selection"] = std::move(selection_obj);
    }
    obj.object["numerics"] = std::move(num_obj);

    // Reference values
    obj.object["ref_length"] = mini_json::value(setup.ref_length);
    obj.object["ref_density"] = mini_json::value(setup.ref_density);
    obj.object["ref_velocity"] = mini_json::value(setup.ref_velocity);
    obj.object["ref_temperature"] = mini_json::value(setup.ref_temperature);
    obj.object["ref_pressure"] = mini_json::value(setup.ref_pressure);

    // Misc
    obj.object["gravity_vector"] = setup.gravity_vector;
    {
        mini_json::array_t bf_arr;
        for (double v : setup.body_force)
            bf_arr.emplace_back(mini_json::value(v));
        obj.object["body_force"] = std::move(bf_arr);
    }
    obj.object["units"] = setup.units;
    obj.object["solver_mode"] = setup.solver_mode;

    mini_json::value sm_obj(mini_json::value_type::object);
    for (const auto& [k, v] : setup.source_metadata) {
        sm_obj.object[k] = v;
    }
    obj.object["source_metadata"] = std::move(sm_obj);
}

// ===========================================================================
// CaseSetup JSON deserialization helpers
// ===========================================================================

static std::string json_get_string(const mini_json::value& obj, const std::string& key,
                                    const std::string& default_val = "") {
    const mini_json::value* v = mini_json::find(obj, key);
    if (v && v->is_string()) return v->string;
    return default_val;
}

static double json_get_number(const mini_json::value& obj, const std::string& key, double default_val = 0.0) {
    const mini_json::value* v = mini_json::find(obj, key);
    if (v && v->is_number()) return v->number;
    return default_val;
}

static bool json_get_bool(const mini_json::value& obj, const std::string& key, bool default_val = false) {
    const mini_json::value* v = mini_json::find(obj, key);
    if (v && v->is_boolean()) return v->boolean;
    return default_val;
}

static BCType bc_type_from_string(const std::string& s) {
    if (s == "wall")            return BCType::WALL;
    if (s == "inlet")           return BCType::INLET;
    if (s == "outlet")          return BCType::OUTLET;
    if (s == "pressure_outlet") return BCType::PRESSURE_OUTLET;
    if (s == "symmetry")        return BCType::SYMMETRY;
    if (s == "periodic")        return BCType::PERIODIC;
    if (s == "interface")       return BCType::INTERFACE;
    if (s == "empty")           return BCType::EMPTY;
    if (s == "internal")        return BCType::INTERNAL;
    return BCType::UNKNOWN;
}

static BCValueType bcvt_from_string(const std::string& s) {
    if (s == "fixed")            return BCValueType::FIXED;
    if (s == "zero_gradient")    return BCValueType::ZERO_GRADIENT;
    if (s == "mixed")            return BCValueType::MIXED;
    if (s == "outlet_pressure")  return BCValueType::OUTLET_PRESSURE;
    if (s == "wall_no_slip")     return BCValueType::WALL_NO_SLIP;
    if (s == "wall_slip")        return BCValueType::WALL_SLIP;
    if (s == "wall_thermal")     return BCValueType::WALL_THERMAL;
    return BCValueType::UNKNOWN;
}

static void case_setup_from_json(const mini_json::value& json, CaseSetup& setup) {
    if (!json.is_object()) return;

    setup.schema_version = static_cast<int>(json_get_number(json, "schema_version", 1.0));

    const mini_json::value* src = mini_json::find(json, "source");
    if (src && src->is_object()) {
        setup.source.solver = json_get_string(*src, "solver");
        setup.source.version = json_get_string(*src, "version");
        setup.source.case_path = json_get_string(*src, "case_path");
        setup.source.case_name = json_get_string(*src, "case_name");
        setup.source.format = json_get_string(*src, "format");
    }

    const mini_json::value* mi = mini_json::find(json, "mesh_info");
    if (mi && mi->is_object()) {
        setup.mesh_info.n_vertices = static_cast<std::size_t>(json_get_number(*mi, "n_vertices"));
        setup.mesh_info.n_faces = static_cast<std::size_t>(json_get_number(*mi, "n_faces"));
        setup.mesh_info.n_cells = static_cast<std::size_t>(json_get_number(*mi, "n_cells"));
        setup.mesh_info.n_boundary_faces = static_cast<std::size_t>(json_get_number(*mi, "n_boundary_faces"));
        setup.mesh_info.n_internal_faces = static_cast<std::size_t>(json_get_number(*mi, "n_internal_faces"));
        setup.mesh_info.n_patches = static_cast<std::size_t>(json_get_number(*mi, "n_patches"));
        setup.mesh_info.dimension = static_cast<int>(json_get_number(*mi, "dimension", 3.0));
        setup.mesh_info.has_polyhedral = json_get_bool(*mi, "has_polyhedral");
        setup.mesh_info.has_nonplanar_faces = json_get_bool(*mi, "has_nonplanar_faces");
        setup.mesh_info.mesh_type = json_get_string(*mi, "mesh_type", "unstructured");

        const mini_json::value* ctypes = mini_json::find(*mi, "cell_types");
        if (ctypes && ctypes->is_array()) {
            for (const auto& ct : ctypes->array) {
                if (ct.is_string()) setup.mesh_info.cell_types.push_back(ct.string);
            }
        }
    }

    setup.physics_model = json_get_string(json, "physics_model");
    setup.turbulence_model = json_get_string(json, "turbulence_model");
    setup.energy_model = json_get_string(json, "energy_model");
    setup.multiphase_model = json_get_string(json, "multiphase_model");
    setup.radiation_model = json_get_string(json, "radiation_model");
    setup.transient = json_get_bool(json, "transient");
    setup.time_step = json_get_number(json, "time_step");
    setup.end_time = json_get_number(json, "end_time");
    setup.max_time_steps = static_cast<int>(json_get_number(json, "max_time_steps"));

    const mini_json::value* mats = mini_json::find(json, "materials");
    if (mats && mats->is_array()) {
        for (const auto& mat_v : mats->array) {
            if (!mat_v.is_object()) continue;
            MaterialSpec mat;
            mat.name = json_get_string(mat_v, "name", "air");
            mat.density = json_get_number(mat_v, "density", 1.225);
            mat.dynamic_viscosity = json_get_number(mat_v, "dynamic_viscosity", 1.789e-5);
            mat.thermal_conductivity = json_get_number(mat_v, "thermal_conductivity", 0.024);
            mat.specific_heat = json_get_number(mat_v, "specific_heat", 1006.43);
            mat.molecular_weight = json_get_number(mat_v, "molecular_weight", 28.97);
            mat.eos_model = json_get_string(mat_v, "eos_model", "ideal_gas");
            mat.turbulence_model = json_get_string(mat_v, "turbulence_model", "laminar");

            const mini_json::value* ep = mini_json::find(mat_v, "extra_properties");
            if (ep && ep->is_object()) {
                for (const auto& [k, v] : ep->object) {
                    if (v.is_string()) mat.extra_properties[k] = v.string;
                    else {
                        std::ostringstream oss;
                        oss << v.number;
                        mat.extra_properties[k] = oss.str();
                    }
                }
            }
            setup.materials.push_back(mat);
        }
    }

    const mini_json::value* bcs = mini_json::find(json, "boundary_conditions");
    if (bcs && bcs->is_array()) {
        for (const auto& bc_v : bcs->array) {
            if (!bc_v.is_object()) continue;
            BoundarySpec bc;
            bc.patch_name = json_get_string(bc_v, "patch_name");
            bc.type = bc_type_from_string(json_get_string(bc_v, "type"));
            bc.value_type = bcvt_from_string(json_get_string(bc_v, "value_type"));
            bc.velocity_magnitude = json_get_number(bc_v, "velocity_magnitude");

            const mini_json::value* vv = mini_json::find(bc_v, "velocity_vector");
            if (vv && vv->is_array()) {
                for (const auto& comp : vv->array) {
                    if (comp.is_number()) bc.velocity_vector.push_back(comp.number);
                }
            }
            bc.temperature = json_get_number(bc_v, "temperature");
            bc.pressure = json_get_number(bc_v, "pressure");
            bc.turbulence_intensity = json_get_number(bc_v, "turbulence_intensity");
            bc.turbulence_length_scale = json_get_number(bc_v, "turbulence_length_scale");
            bc.roughness_height = json_get_number(bc_v, "roughness_height");
            bc.wall_temperature = json_get_number(bc_v, "wall_temperature");
            bc.back_pressure = json_get_number(bc_v, "back_pressure");
            bc.source_zone_name = json_get_string(bc_v, "source_zone_name");
            bc.source_type_name = json_get_string(bc_v, "source_type_name");

            const mini_json::value* rp = mini_json::find(bc_v, "raw_params");
            if (rp && rp->is_object()) {
                for (const auto& [k, v] : rp->object) {
                    if (v.is_string()) bc.raw_params[k] = v.string;
                    else {
                        std::ostringstream oss;
                        oss << v.number;
                        bc.raw_params[k] = oss.str();
                    }
                }
            }
            setup.boundary_conditions.push_back(bc);
        }
    }

    const mini_json::value* ic = mini_json::find(json, "initial_condition");
    if (ic && ic->is_object()) {
        setup.has_initial_condition = true;
        setup.initial_condition.velocity = json_get_number(*ic, "velocity");
        setup.initial_condition.pressure = json_get_number(*ic, "pressure");
        setup.initial_condition.temperature = json_get_number(*ic, "temperature", 288.15);

        const mini_json::value* vv = mini_json::find(*ic, "velocity_vector");
        if (vv && vv->is_array()) {
            for (const auto& comp : vv->array) {
                if (comp.is_number()) setup.initial_condition.velocity_vector.push_back(comp.number);
            }
        }

        const mini_json::value* sf = mini_json::find(*ic, "scalar_fields");
        if (sf && sf->is_object()) {
            for (const auto& [k, v] : sf->object) {
                if (v.is_number()) setup.initial_condition.scalar_fields[k] = v.number;
            }
        }
    }

    const mini_json::value* num = mini_json::find(json, "numerics");
    if (num && num->is_object()) {
        setup.numerics.momentum_scheme = json_get_string(*num, "momentum_scheme", "first_order");
        setup.numerics.pressure_scheme = json_get_string(*num, "pressure_scheme", "standard");
        setup.numerics.momentum_interpolation = json_get_string(*num, "momentum_interpolation", "linear");
        setup.numerics.transient_scheme = json_get_string(*num, "transient_scheme", "steady");
        setup.numerics.gradient_operator = json_get_string(*num, "gradient_operator", "green_gauss_cell");
        setup.numerics.under_relaxation_momentum = json_get_number(*num, "under_relaxation_momentum", 0.7);
        setup.numerics.under_relaxation_pressure = json_get_number(*num, "under_relaxation_pressure", 0.3);
        setup.numerics.coupled_solver = json_get_string(*num, "coupled_solver", "SIMPLE");
        setup.numerics.preconditioner = json_get_string(*num, "preconditioner");
        setup.numerics.residual_target = json_get_string(*num, "residual_target", "1e-5");
        setup.numerics.max_iterations = static_cast<int>(json_get_number(*num, "max_iterations", 500.0));

        const mini_json::value* selection = mini_json::find(*num, "selection");
        if (selection && selection->is_object()) {
            setup.has_explicit_numerics = true;
            const mini_json::value* entries = mini_json::find(*selection, "entries");
            if (entries && entries->is_array()) {
                for (const auto& entry_v : entries->array) {
                    if (!entry_v.is_object()) {
                        setup.numerical_report.errors.push_back(
                            "numerical selection entry must be an object");
                        continue;
                    }
                    const std::string family = json_get_string(entry_v, "family");
                    const std::string key = json_get_string(entry_v, "configuration_key");
                    cfdx::core::NumericalMethodFamily parsed = cfdx::core::NumericalMethodFamily::Convection;
                    bool known_family = false;
                    for (const auto candidate : {
                        cfdx::core::NumericalMethodFamily::Gradient,
                        cfdx::core::NumericalMethodFamily::Interpolation,
                        cfdx::core::NumericalMethodFamily::Convection,
                        cfdx::core::NumericalMethodFamily::Diffusion,
                        cfdx::core::NumericalMethodFamily::Source,
                        cfdx::core::NumericalMethodFamily::Temporal,
                        cfdx::core::NumericalMethodFamily::TimeStep,
                        cfdx::core::NumericalMethodFamily::Nonlinear,
                        cfdx::core::NumericalMethodFamily::LinearSolver,
                        cfdx::core::NumericalMethodFamily::Preconditioner,
                        cfdx::core::NumericalMethodFamily::PressureVelocity,
                        cfdx::core::NumericalMethodFamily::Conservation,
                        cfdx::core::NumericalMethodFamily::Reconstruction,
                        cfdx::core::NumericalMethodFamily::Schur}) {
                        if (family == cfdx::core::to_string(candidate)) {
                            parsed = candidate;
                            known_family = true;
                            break;
                        }
                    }
                    if (known_family)
                        setup.numerical_config.entries.push_back({parsed, key});
                    else
                        setup.numerical_report.errors.push_back(
                            "unknown numerical method family: " + family);
                }
            }
            const mini_json::value* required = mini_json::find(*selection, "required_families");
            if (required && required->is_array()) {
                for (const auto& family_v : required->array) {
                    if (!family_v.is_string()) {
                        setup.numerical_report.errors.push_back(
                            "required numerical method family must be a string");
                        continue;
                    }
                    bool known_family = false;
                    for (const auto candidate : {
                        cfdx::core::NumericalMethodFamily::Gradient,
                        cfdx::core::NumericalMethodFamily::Interpolation,
                        cfdx::core::NumericalMethodFamily::Convection,
                        cfdx::core::NumericalMethodFamily::Diffusion,
                        cfdx::core::NumericalMethodFamily::Source,
                        cfdx::core::NumericalMethodFamily::Temporal,
                        cfdx::core::NumericalMethodFamily::TimeStep,
                        cfdx::core::NumericalMethodFamily::Nonlinear,
                        cfdx::core::NumericalMethodFamily::LinearSolver,
                        cfdx::core::NumericalMethodFamily::Preconditioner,
                        cfdx::core::NumericalMethodFamily::PressureVelocity,
                        cfdx::core::NumericalMethodFamily::Conservation,
                        cfdx::core::NumericalMethodFamily::Reconstruction,
                        cfdx::core::NumericalMethodFamily::Schur}) {
                        if (family_v.string == cfdx::core::to_string(candidate)) {
                            setup.numerical_config.required_families.push_back(candidate);
                            known_family = true;
                            break;
                        }
                    }
                    if (!known_family)
                        setup.numerical_report.errors.push_back(
                            "unknown required numerical method family: " + family_v.string);
                }
            }
        }

        const mini_json::value* rs = mini_json::find(*num, "raw_settings");
        if (rs && rs->is_object()) {
            for (const auto& [k, v] : rs->object) {
                if (v.is_string()) setup.numerics.raw_settings[k] = v.string;
                else {
                    std::ostringstream oss;
                    oss << v.number;
                    setup.numerics.raw_settings[k] = oss.str();
                }
            }
        }
    }

    setup.ref_length = json_get_number(json, "ref_length", 1.0);
    setup.ref_density = json_get_number(json, "ref_density", 1.0);
    setup.ref_velocity = json_get_number(json, "ref_velocity", 1.0);
    setup.ref_temperature = json_get_number(json, "ref_temperature", 293.15);
    setup.ref_pressure = json_get_number(json, "ref_pressure", 101325.0);

    setup.gravity_vector = json_get_string(json, "gravity_vector");
    const mini_json::value* bf = mini_json::find(json, "body_force");
    if (bf && bf->is_array()) {
        for (const auto& comp : bf->array) {
            if (comp.is_number()) setup.body_force.push_back(comp.number);
        }
    }
    setup.units = json_get_string(json, "units", "SI");
    setup.solver_mode = json_get_string(json, "solver_mode", "steady");

    const mini_json::value* sm = mini_json::find(json, "source_metadata");
    if (sm && sm->is_object()) {
        for (const auto& [k, v] : sm->object) {
            if (v.is_string()) setup.source_metadata[k] = v.string;
            else {
                std::ostringstream oss;
                oss << v.number;
                setup.source_metadata[k] = oss.str();
            }
        }
    }
}

// ===========================================================================
// GapAnalysis JSON deserialization
// ===========================================================================

static Severity severity_from_string(const std::string& s) {
    if (s == "supported")                  return Severity::SUPPORTED;
    if (s == "approximated")               return Severity::APPROXIMATED;
    if (s == "unsupported_nonblocking")    return Severity::UNSUPPORTED_NONBLOCK;
    if (s == "unsupported_blocking")       return Severity::UNSUPPORTED_BLOCK;
    if (s == "unavailable")                return Severity::UNAVAILABLE;
    return Severity::SUPPORTED;
}

static void gap_analysis_from_json(const mini_json::value& json, GapAnalysis& gap) {
    if (!json.is_object()) return;

    const mini_json::value* findings = mini_json::find(json, "findings");
    if (findings && findings->is_array()) {
        for (const auto& f : findings->array) {
            if (!f.is_object()) continue;
            Finding finding;
            finding.severity = severity_from_string(json_get_string(f, "severity"));
            finding.category = json_get_string(f, "category");
            finding.feature = json_get_string(f, "feature");
            finding.detail = json_get_string(f, "detail");
            finding.suggestion = json_get_string(f, "suggestion");
            gap.add(finding.severity, finding.category, finding.feature,
                    finding.detail, finding.suggestion);
        }
    }
}

// ===========================================================================
// read_case_cfdx_h5
// ===========================================================================

bool read_case_cfdx_h5(const std::string& filename,
                       cfdx::core::Mesh& mesh,
                       SourceInfo& source,
                       CaseSetup& setup,
                       GapAnalysis& gap) {
    hid_t fapl = case_create_file_access_plist();
    if (fapl < 0) return false;
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);
    H5Pclose(fapl);
    if (file < 0) return false;

    auto fail = [&](const std::string& reason) {
        std::cerr << "HDF5 case IO error in '" << filename
                  << "': " << reason << '\n';
        H5Fclose(file);
        return false;
    };

    // --- Source info (attributes that may be absent in mesh-only files) ---
    case_read_attr_str(file, "source_solver", source.solver);
    case_read_attr_str(file, "source_format", source.format);
    case_read_attr_str(file, "source_version", source.version);
    case_read_attr_str(file, "source_case_path", source.case_path);
    case_read_attr_str(file, "source_case_name", source.case_name);

    // --- Case setup JSON (if present) ---
    std::string case_setup_json;
    if (case_read_attr_str(file, "case_setup_json", case_setup_json)) {
        mini_json::value root = mini_json::value::parse_safe(case_setup_json);
        if (root.is_object()) {
            case_setup_from_json(root, setup);
        }
    }

    // Resolve the canonical numerical selection as part of case loading.
    // A malformed explicit selection is a case-load failure, not a solver-time
    // fallback. Legacy files without the block remain readable during migration.
    if (setup.has_explicit_numerics) {
        const auto parse_errors = setup.numerical_report.errors;
        setup.numerical_report = cfdx::core::resolve_case_numerics(
            source.case_name.empty() ? filename : source.case_name,
            setup.numerical_config);
        setup.numerical_report.errors.insert(
            setup.numerical_report.errors.end(), parse_errors.begin(), parse_errors.end());
        if (!setup.numerical_report.valid()) {
            std::string reason = "invalid explicit numerical configuration";
            if (!setup.numerical_report.errors.empty())
                reason += ": " + setup.numerical_report.errors.front();
            return fail(reason);
        }
    }

    // --- Gap report JSON (if present) ---
    std::string gap_json;
    if (case_read_attr_str(file, "gap_report_json", gap_json)) {
        mini_json::value root = mini_json::value::parse_safe(gap_json);
        if (root.is_object()) {
            gap_analysis_from_json(root, gap);
        }
    }

    // --- Mesh topology ---
    // Read integrity attributes (may be absent in Python-written files)
    std::string format_version, schema_version, topology_hash, geometry_hash, mesh_hash;
    const bool has_format = case_read_attr_str(file, "format_version", format_version);
    const bool has_schema = case_read_attr_str(file, "schema_version", schema_version);
    const bool has_topology = case_read_attr_str(file, "topology_hash", topology_hash);
    const bool has_geometry = case_read_attr_str(file, "geometry_hash", geometry_hash);
    const bool has_mesh = case_read_attr_str(file, "mesh_hash", mesh_hash);

    // Validate versions if present
    if (has_format || has_schema) {
        if (has_format) {
            try {
                if (std::stoul(format_version) != CFDX_HDF5_FORMAT_VERSION)
                    return fail("unsupported HDF5 format version " + format_version);
            } catch (...) {
                return fail("invalid HDF5 format version");
            }
        }
        if (has_schema) {
            try {
                if (std::stoul(schema_version) != CFDX_HDF5_SCHEMA_VERSION)
                    return fail("unsupported HDF5 schema version " + schema_version);
            } catch (...) {
                return fail("invalid HDF5 schema version");
            }
        }
    }

    // Read core topology datasets
    std::vector<double> pts;
    std::vector<std::uint64_t> fv, fo, owner, cf, co;
    std::vector<std::int64_t> neighbour;

    if (!case_read_dataset_double(file, "points", pts, 2) ||
        !case_read_dataset_u64(file, "face_vertices", fv) ||
        !case_read_dataset_u64(file, "face_offsets", fo) ||
        !case_read_dataset_u64(file, "owner", owner) ||
        !case_read_dataset_i64(file, "neighbour", neighbour) ||
        !case_read_dataset_u64(file, "cell_faces", cf) ||
        !case_read_dataset_u64(file, "cell_offsets", co)) {
        return fail("missing required topology dataset");
    }

    if (pts.size() % 3 != 0 || fo.empty() || co.empty() ||
        fo.front() != 0 || co.front() != 0 ||
        fo.back() != fv.size() || co.back() != cf.size() ||
        owner.size() != neighbour.size() ||
        fo.size() - 1 != owner.size()) {
        return fail("invalid topology dataset dimensions or CSR terminal offsets");
    }

    // Validate CSR offsets
    for (std::size_t i = 1; i < fo.size(); ++i) {
        if (fo[i] < fo[i - 1] || fo[i] > fv.size())
            return fail("face CSR offsets are not monotonic or exceed face-vertex storage");
    }
    for (std::size_t i = 1; i < co.size(); ++i) {
        if (co[i] < co[i - 1] || co[i] > cf.size())
            return fail("cell CSR offsets are not monotonic or exceed cell-face storage");
    }

    const std::size_t n_points = pts.size() / 3;
    const std::size_t n_faces = fo.size() - 1;
    const std::size_t n_cells = co.size() - 1;

    // Validate integrity hashes ONLY when ALL three are present.
    // Python-written files (mesh_topology attribute but no hashes) are accepted.
    const bool has_all_hashes = has_topology && has_mesh && has_geometry;
    const bool has_any_hash = has_topology || has_geometry || has_mesh;

    if (has_any_hash && !has_all_hashes) {
        return fail("incomplete integrity hashes (some present, some absent)");
    }

    if (has_all_hashes) {
        std::uint64_t topology = 1469598103934665603ULL;
        topology = case_fnv1a_update_vector(topology, fv.data(), fv.size());
        if (n_faces > 0) topology = case_fnv1a_update_vector(topology, fo.data(), fo.size());
        topology = case_fnv1a_update_vector(topology, owner.data(), owner.size());
        topology = case_fnv1a_update_vector(topology, neighbour.data(), neighbour.size());
        topology = case_fnv1a_update_vector(topology, cf.data(), cf.size());
        if (n_cells > 0) topology = case_fnv1a_update_vector(topology, co.data(), co.size());
        if (topology_hash != case_hash_hex(topology))
            return fail("topology integrity hash mismatch");

        std::vector<double> xs, ys, zs;
        xs.reserve(n_points); ys.reserve(n_points); zs.reserve(n_points);
        for (std::size_t i = 0; i < n_points; ++i) {
            xs.push_back(pts[i * 3]);
            ys.push_back(pts[i * 3 + 1]);
            zs.push_back(pts[i * 3 + 2]);
        }

        std::uint64_t geometry = 1469598103934665603ULL;
        geometry = case_fnv1a_update_vector(geometry, xs.data(), xs.size());
        geometry = case_fnv1a_update_vector(geometry, ys.data(), ys.size());
        geometry = case_fnv1a_update_vector(geometry, zs.data(), zs.size());
        if (geometry_hash != case_hash_hex(geometry))
            return fail("geometry integrity hash mismatch");

        std::uint64_t full_mesh = topology;
        full_mesh = case_fnv1a_update_vector(full_mesh, xs.data(), xs.size());
        full_mesh = case_fnv1a_update_vector(full_mesh, ys.data(), ys.size());
        full_mesh = case_fnv1a_update_vector(full_mesh, zs.data(), zs.size());
        if (mesh_hash != case_hash_hex(full_mesh))
            return fail("mesh integrity hash mismatch");
    }

    // Validate vertex indices
    for (std::size_t j = 0; j < fv.size(); ++j) {
        if (fv[j] >= n_points)
            return fail("face vertex index is outside the point array");
    }

    // Build mesh
    mesh.clear();
    mesh.points().resize(n_points);
    for (std::size_t i = 0; i < n_points; ++i) {
        mesh.points().set(i, pts[i * 3], pts[i * 3 + 1], pts[i * 3 + 2]);
    }

    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::uint64_t begin = fo[f];
        const std::uint64_t end = fo[f + 1];
        std::vector<cfdx::core::FaceConnectivity::Index> vertices;
        vertices.reserve(static_cast<std::size_t>(end - begin));
        for (std::uint64_t j = begin; j < end; ++j) {
            vertices.push_back(static_cast<cfdx::core::FaceConnectivity::Index>(fv[j]));
        }
        if (vertices.size() < 3)
            return fail("face contains fewer than three vertices");
        mesh.faces().push_face(std::move(vertices));
    }

    mesh.ownership().resize(n_faces);
    for (std::size_t i = 0; i < n_faces; ++i) {
        if (owner[i] >= n_cells)
            return fail("owner index is outside the cell range");
        if (neighbour[i] < -1 ||
            (neighbour[i] >= 0 && static_cast<std::uint64_t>(neighbour[i]) >= n_cells)) {
            return fail("neighbour index is outside the cell range");
        }
        mesh.ownership().set_owner(i, owner[i]);
        mesh.ownership().set_neighbour(i, neighbour[i]);
    }

    for (std::size_t c = 0; c < n_cells; ++c) {
        const std::uint64_t begin = co[c];
        const std::uint64_t end = co[c + 1];
        std::vector<cfdx::core::CellConnectivity::FaceId> faces;
        faces.reserve(static_cast<std::size_t>(end - begin));
        for (std::uint64_t j = begin; j < end; ++j) {
            if (cf[j] >= n_faces)
                return fail("cell-face index is outside the face range");
            faces.push_back(static_cast<cfdx::core::CellConnectivity::FaceId>(cf[j]));
        }
        if (faces.empty())
            return fail("cell contains no faces");
        mesh.cells().push_cell(std::move(faces));
    }

    // Boundary patches
    std::string patches_str;
    if (case_read_attr_str(file, "boundary_patches", patches_str)) {
        const auto patch_entries = case_parse_patch_metadata(patches_str);
        std::vector<std::uint64_t> patch_face_ids, patch_face_offsets;
        if (!case_read_dataset_u64(file, "patch_face_ids", patch_face_ids) ||
            !case_read_dataset_u64(file, "patch_face_offsets", patch_face_offsets) ||
            patch_face_offsets.empty() ||
            patch_face_offsets.front() != 0 ||
            patch_face_offsets.back() != patch_face_ids.size() ||
            patch_face_offsets.size() != patch_entries.size() + 1) {
            return fail("boundary patch metadata datasets are missing or inconsistent");
        }

        for (std::size_t i = 1; i < patch_face_offsets.size(); ++i) {
            if (patch_face_offsets[i] < patch_face_offsets[i - 1] ||
                patch_face_offsets[i] > patch_face_ids.size()) {
                return fail("boundary patch CSR offsets are invalid");
            }
        }

        cfdx::core::BoundaryPatches bp;
        for (std::size_t i = 0; i < patch_entries.size(); ++i) {
            const std::string& entry = patch_entries[i];
            const size_t p1 = entry.find(':');
            const size_t p2 = entry.find(':', p1 == std::string::npos ? 0 : p1 + 1);
            const size_t p3 = entry.find(':', p2 == std::string::npos ? 0 : p2 + 1);
            if (p1 == std::string::npos || p2 == std::string::npos || p3 == std::string::npos) {
                return fail("boundary patch metadata entry is malformed");
            }

            cfdx::core::Patch patch;
            patch.name = entry.substr(0, p1);
            try {
                const std::uint64_t start = std::stoull(entry.substr(p1 + 1, p2 - p1 - 1));
                const std::uint64_t count = std::stoull(entry.substr(p2 + 1, p3 - p2 - 1));
                patch.type = static_cast<cfdx::core::PatchType>(
                    std::stoull(entry.substr(p3 + 1)));

                const std::uint64_t offset_begin = patch_face_offsets[i];
                const std::uint64_t offset_end = patch_face_offsets[i + 1];
                if (count != offset_end - offset_begin ||
                    start > n_faces || start + count > n_faces) {
                    return fail("boundary patch face count/range is inconsistent");
                }

                patch.face_ids.reserve(static_cast<std::size_t>(count));
                for (std::uint64_t j = 0; j < count; ++j) {
                    const std::uint64_t face_id = patch_face_ids[offset_begin + j];
                    if (face_id >= n_faces)
                        return fail("boundary patch references a face outside the mesh");
                    patch.face_ids.push_back(static_cast<cfdx::core::FaceIndex>(face_id));
                }
            } catch (...) {
                return fail("mesh topology validation failed");
            }
            bp.add_patch(std::move(patch));
        }
        mesh.set_boundary(bp);
    } else {
        // Reconstruct conservative generic boundary from face ownership.
        cfdx::core::BoundaryPatches bp;
        cfdx::core::Patch patch;
        patch.name = "boundary";
        patch.type = cfdx::core::PatchType::WALL;
        for (std::size_t f = 0; f < neighbour.size(); ++f) {
            if (neighbour[f] < 0)
                patch.face_ids.push_back(static_cast<cfdx::core::FaceIndex>(f));
        }
        if (!patch.face_ids.empty())
            bp.add_patch(std::move(patch));
        mesh.set_boundary(bp);
    }

    const bool valid = mesh.topo_validate().ok;
    H5Fclose(file);
    return valid;
}

// ===========================================================================
// write_case_cfdx_h5
// ===========================================================================

static std::string utc_timestamp() {
    const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &now);
#else
    gmtime_r(&now, &tm);
#endif
    std::ostringstream os;
    os << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
    return os.str();
}

static void write_schema_attributes(hid_t file, const cfdx::core::Mesh& mesh) {
    case_write_attr_str(file, "format_version", std::to_string(CFDX_HDF5_FORMAT_VERSION));
    case_write_attr_str(file, "schema_version", std::to_string(CFDX_HDF5_SCHEMA_VERSION));
    case_write_attr_str(file, "cfdx_version", CFDX_VERSION);

    std::uint64_t topology = 1469598103934665603ULL;
    topology = case_fnv1a_update_vector(topology, mesh.faces().vertices_data(), mesh.faces().n_vertices());
    topology = case_fnv1a_update_vector(topology, mesh.faces().offsets_data(), mesh.faces().n_faces() + 1);
    topology = case_fnv1a_update_vector(topology, mesh.ownership().owner_data(), mesh.ownership().size());
    topology = case_fnv1a_update_vector(topology, mesh.ownership().neighbour_data(), mesh.ownership().size());
    topology = case_fnv1a_update_vector(topology, mesh.cells().faces_data(), mesh.cells().n_face_refs());
    topology = case_fnv1a_update_vector(topology, mesh.cells().offsets_data(), mesh.cells().n_cells() + 1);
    case_write_attr_str(file, "topology_hash", case_hash_hex(topology));

    std::uint64_t geometry = 1469598103934665603ULL;
    geometry = case_fnv1a_update_vector(geometry, mesh.points().x_data(), mesh.n_points());
    geometry = case_fnv1a_update_vector(geometry, mesh.points().y_data(), mesh.n_points());
    geometry = case_fnv1a_update_vector(geometry, mesh.points().z_data(), mesh.n_points());
    case_write_attr_str(file, "geometry_hash", case_hash_hex(geometry));

    std::uint64_t full_mesh = topology;
    full_mesh = case_fnv1a_update_vector(full_mesh, mesh.points().x_data(), mesh.n_points());
    full_mesh = case_fnv1a_update_vector(full_mesh, mesh.points().y_data(), mesh.n_points());
    full_mesh = case_fnv1a_update_vector(full_mesh, mesh.points().z_data(), mesh.n_points());
    case_write_attr_str(file, "mesh_hash", case_hash_hex(full_mesh));

    case_write_attr_str(file, "creation_date", utc_timestamp());
    case_write_attr_str(file, "modification_date", utc_timestamp());
    case_write_attr_str(file, "dimension", "3");
    case_write_attr_str(file, "precision", "float64");
    case_write_attr_str(file, "endian", "native");
    case_write_attr_str(file, "mesh_topology", "cfdx-csr-v1");
}

bool write_case_cfdx_h5(const std::string& filename,
                        const cfdx::core::Mesh& mesh,
                        const SourceInfo& source,
                        const CaseSetup& setup,
                        const GapAnalysis& gap) {
    cfdx::core::CaseNumericsReport numerical_report = setup.numerical_report;
    if (setup.has_explicit_numerics) {
        numerical_report = cfdx::core::resolve_case_numerics(
            source.case_name.empty() ? filename : source.case_name,
            setup.numerical_config);
        if (!numerical_report.valid()) return false;
    }

    hid_t fapl = case_create_file_access_plist();
    if (fapl < 0) return false;
    hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    H5Pclose(fapl);
    if (file < 0) return false;

    // --- Schema + integrity hashes ---
    write_schema_attributes(file, mesh);

    // --- Mesh counts ---
    case_write_attr_str(file, "n_points", std::to_string(mesh.n_points()));
    case_write_attr_str(file, "n_faces", std::to_string(mesh.n_faces()));
    case_write_attr_str(file, "n_cells", std::to_string(mesh.n_cells()));

    // --- Source info ---
    case_write_attr_str(file, "source_solver", source.solver);
    case_write_attr_str(file, "source_format", source.format);
    case_write_attr_str(file, "source_version", source.version);
    case_write_attr_str(file, "source_case_path", source.case_path);
    case_write_attr_str(file, "source_case_name", source.case_name);

    // --- Case setup JSON ---
    mini_json::value setup_json;
    case_setup_to_json(setup, setup_json);
    case_write_attr_str(file, "case_setup_json", setup_json.serialize());

    // --- Gap report JSON ---
    case_write_attr_str(file, "gap_report_json", gap.to_json());
    if (setup.has_explicit_numerics) {
        case_write_attr_str(file, "numerical_selection_report",
                       cfdx::core::format_numerics_report(numerical_report));
    }

    // --- Points ---
    {
        const hsize_t dims[2] = {mesh.n_points(), 3};
        std::vector<double> pts(mesh.n_points() * 3);
        const cfdx::core::PointCloud& pc = mesh.points();
        for (std::size_t i = 0; i < mesh.n_points(); ++i) {
            pts[i * 3 + 0] = pc.x(i);
            pts[i * 3 + 1] = pc.y(i);
            pts[i * 3 + 2] = pc.z(i);
        }
        if (dims[0] > 0) {
            write_dataset_double_2d(file, "points", pts.data(), dims[0], 3);
        } else {
            const hsize_t zero[2] = {0, 3};
            hid_t space = H5Screate_simple(2, zero, nullptr);
            if (space >= 0) {
                hid_t ds = H5Dcreate2(file, "points", H5T_NATIVE_DOUBLE, space,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (ds >= 0) H5Dclose(ds);
                H5Sclose(space);
            }
        }
    }

    // --- Faces (CSR) ---
    {
        const cfdx::core::FaceConnectivity& fc = mesh.faces();
        case_write_dataset_u64(file, "face_vertices", fc.vertices_data(), fc.n_vertices());
        case_write_dataset_u64(file, "face_offsets", fc.offsets_data(), fc.n_faces() + 1);
    }

    // --- Owner / neighbour ---
    {
        const cfdx::core::FaceOwnership& own = mesh.ownership();
        case_write_dataset_u64(file, "owner", own.owner_data(), own.size());
        case_write_dataset_i64(file, "neighbour", own.neighbour_data(), own.size());
    }

    // --- Cell faces (CSR) ---
    {
        const cfdx::core::CellConnectivity& cc = mesh.cells();
        case_write_dataset_u64(file, "cell_faces", cc.faces_data(), cc.n_face_refs());
        case_write_dataset_u64(file, "cell_offsets", cc.offsets_data(), cc.n_cells() + 1);
    }

    // --- Boundary patches ---
    {
        const cfdx::core::BoundaryPatches& bp = mesh.boundary();
        std::size_t n_patches = bp.n_patches();
        if (n_patches > 0) {
            std::string patches_str;
            for (std::size_t i = 0; i < n_patches; ++i) {
                const auto& patch = bp.patch(i);
                if (i > 0) patches_str += ";";
                patches_str += patch.name + ":" +
                    std::to_string(patch.face_ids.empty() ? 0 : patch.face_ids[0]) + ":" +
                    std::to_string(patch.face_ids.size()) + ":" +
                    std::to_string(static_cast<std::uint32_t>(patch.type));
            }
            case_write_attr_str(file, "boundary_patches", patches_str);

            std::vector<std::uint64_t> all_face_ids;
            for (std::size_t i = 0; i < n_patches; ++i) {
                const auto& patch = bp.patch(i);
                all_face_ids.insert(all_face_ids.end(), patch.face_ids.begin(), patch.face_ids.end());
            }
            std::vector<std::uint64_t> patch_face_offsets(n_patches + 1);
            patch_face_offsets[0] = 0;
            for (std::size_t i = 0; i < n_patches; ++i) {
                patch_face_offsets[i + 1] = patch_face_offsets[i] + bp.patch(i).face_ids.size();
            }
            case_write_dataset_u64(file, "patch_face_ids", all_face_ids.data(), all_face_ids.size());
            case_write_dataset_u64(file, "patch_face_offsets", patch_face_offsets.data(), patch_face_offsets.size());
        }
    }

    H5Fclose(file);
    return true;
}

// ===========================================================================
// Field read/write helpers
// ===========================================================================

static bool case_ensure_group(hid_t file_id, const std::string& path) {
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

bool read_scalar_fields_hdf5(const std::string& filename,
                             std::vector<std::pair<std::string, cfdx::core::ScalarCellField>>& fields) {
    hid_t fapl = case_create_file_access_plist();
    if (fapl < 0) return false;
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);
    H5Pclose(fapl);
    if (file < 0) return false;

    hid_t fields_grp = H5Gopen2(file, "fields", H5P_DEFAULT);
    if (fields_grp < 0) { H5Fclose(file); return false; }

    hid_t scalar_grp = H5Gopen2(fields_grp, "scalar", H5P_DEFAULT);
    if (scalar_grp < 0) { H5Gclose(fields_grp); H5Fclose(file); return true; }

    hsize_t num_objs = 0;
    H5Gget_num_objs(scalar_grp, &num_objs);

    bool success = true;
    for (hsize_t i = 0; i < num_objs; ++i) {
        char obj_name[256];
        H5Gget_objname_by_idx(scalar_grp, i, obj_name, sizeof(obj_name));
        std::string name(obj_name);

        hid_t ds = H5Dopen2(scalar_grp, name.c_str(), H5P_DEFAULT);
        if (ds < 0) { success = false; continue; }

        hid_t space = H5Dget_space(ds);
        hsize_t dims[1] = {0};
        H5Sget_simple_extent_dims(space, dims, nullptr);

        std::vector<double> values(static_cast<std::size_t>(dims[0]));
        if (dims[0] > 0) {
            H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data());
        }

        cfdx::core::ScalarCellField field(static_cast<std::size_t>(dims[0]), name, "unit", 1);
        for (std::size_t j = 0; j < values.size(); ++j) {
            field(static_cast<std::size_t>(j)) = values[j];
        }
        fields.emplace_back(name, std::move(field));

        H5Sclose(space);
        H5Dclose(ds);
    }

    H5Gclose(scalar_grp);
    H5Gclose(fields_grp);
    H5Fclose(file);
    return success;
}

bool read_vector_fields_hdf5(const std::string& filename,
                             std::vector<std::pair<std::string, cfdx::core::Vec3CellField>>& fields) {
    hid_t fapl = case_create_file_access_plist();
    if (fapl < 0) return false;
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);
    H5Pclose(fapl);
    if (file < 0) return false;

    hid_t fields_grp = H5Gopen2(file, "fields", H5P_DEFAULT);
    if (fields_grp < 0) { H5Fclose(file); return false; }

    hid_t vec_grp = H5Gopen2(fields_grp, "vector", H5P_DEFAULT);
    if (vec_grp < 0) { H5Gclose(fields_grp); H5Fclose(file); return true; }

    hsize_t num_objs = 0;
    H5Gget_num_objs(vec_grp, &num_objs);

    bool success = true;
    for (hsize_t i = 0; i < num_objs; ++i) {
        char obj_name[256];
        H5Gget_objname_by_idx(vec_grp, i, obj_name, sizeof(obj_name));
        std::string name(obj_name);

        hid_t ds = H5Dopen2(vec_grp, name.c_str(), H5P_DEFAULT);
        if (ds < 0) { success = false; continue; }

        hid_t space = H5Dget_space(ds);
        int rank = H5Sget_simple_extent_ndims(space);
        hsize_t dims[2] = {0, 0};
        H5Sget_simple_extent_dims(space, dims, nullptr);

        std::size_t n_cells = static_cast<std::size_t>(dims[0]);
        std::size_t dim = rank == 2 ? static_cast<std::size_t>(dims[1]) : 1;

        cfdx::core::Vec3CellField field(n_cells, name, "m/s", 3);

        if (n_cells > 0 && dim >= 3) {
            std::vector<double> flat(n_cells * dim);
            H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat.data());
            for (std::size_t c = 0; c < n_cells; ++c) {
                field.set(c, flat[c * dim], flat[c * dim + 1], flat[c * dim + 2]);
            }
        } else if (n_cells > 0) {
            std::vector<double> flat(n_cells * dim);
            H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat.data());
            for (std::size_t c = 0; c < n_cells; ++c) {
                field.set(c, flat[c * dim], 0.0, 0.0);
            }
        }

        fields.emplace_back(name, std::move(field));

        H5Sclose(space);
        H5Dclose(ds);
    }

    H5Gclose(vec_grp);
    H5Gclose(fields_grp);
    H5Fclose(file);
    return success;
}

bool write_fields_hdf5(const std::string& filename,
                       const std::vector<cfdx::core::ScalarCellField>& scalars,
                       const std::vector<cfdx::core::Vec3CellField>& vectors) {
    hid_t fapl = case_create_file_access_plist();
    if (fapl < 0) return false;
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl);
    H5Pclose(fapl);
    if (file < 0) return false;

    if (!case_ensure_group(file, "fields/scalar")) {
        H5Fclose(file);
        return false;
    }
    if (!case_ensure_group(file, "fields/vector")) {
        H5Fclose(file);
        return false;
    }

    hid_t scalar_grp = H5Gopen2(file, "fields/scalar", H5P_DEFAULT);
    if (scalar_grp < 0) { H5Fclose(file); return false; }

    for (const auto& sf : scalars) {
        if (sf.empty()) continue;
        hsize_t n = sf.size();
        hid_t space = H5Screate_simple(1, &n, nullptr);
        if (space < 0) continue;
        hid_t ds = H5Dcreate2(scalar_grp, sf.name().c_str(), H5T_NATIVE_DOUBLE, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (ds >= 0) {
            std::vector<double> flat(n);
            for (std::size_t i = 0; i < n; ++i) {
                flat[i] = sf(static_cast<std::size_t>(i));
            }
            H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat.data());
            H5Dclose(ds);
        }
        H5Sclose(space);
    }
    H5Gclose(scalar_grp);

    hid_t vec_grp = H5Gopen2(file, "fields/vector", H5P_DEFAULT);
    if (vec_grp < 0) { H5Fclose(file); return false; }

    for (const auto& vf : vectors) {
        if (vf.empty()) continue;
        hsize_t dims[2] = {vf.size(), 3};
        hid_t space = H5Screate_simple(2, dims, nullptr);
        if (space < 0) continue;
        hid_t ds = H5Dcreate2(vec_grp, vf.name().c_str(), H5T_NATIVE_DOUBLE, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (ds >= 0) {
            std::vector<double> flat(vf.size() * 3);
            for (std::size_t i = 0; i < vf.size(); ++i) {
                flat[i * 3 + 0] = vf(i, 0);
                flat[i * 3 + 1] = vf(i, 1);
                flat[i * 3 + 2] = vf(i, 2);
            }
            H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, flat.data());
            H5Dclose(ds);
        }
        H5Sclose(space);
    }
    H5Gclose(vec_grp);

    H5Fclose(file);
    return true;
}

}  // namespace io
}  // namespace cfdx

