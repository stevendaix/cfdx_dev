// M0.9-T02 — HDF5 reader

#include "hdf5_reader.h"

#include <H5public.h>
#include <H5Dpublic.h>
#include <H5Fpublic.h>
#include <H5Gpublic.h>
#include <H5Tpublic.h>
#include <H5Spublic.h>
#include <H5Apublic.h>
#include <H5Ppublic.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

namespace cfdx {
namespace io {

static bool read_dataset_double(hid_t loc_id, const char* name,
                                std::vector<double>& out) {
    hid_t ds = H5Dopen2(loc_id, name, H5P_DEFAULT);
    if (ds < 0) return false;
    hid_t space = H5Dget_space(ds);
    hsize_t dims[2];
    int rank = H5Sget_simple_extent_dims(space, dims, nullptr);
    hsize_t total = 1;
    for (int i = 0; i < rank; ++i) total *= dims[i];
    out.resize(total);
    herr_t status = H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
    H5Sclose(space);
    H5Dclose(ds);
    return status >= 0;
}

static bool read_dataset_u64(hid_t loc_id, const char* name,
                             std::vector<std::uint64_t>& out) {
    hid_t ds = H5Dopen2(loc_id, name, H5P_DEFAULT);
    if (ds < 0) return false;
    hid_t space = H5Dget_space(ds);
    hsize_t dims[2];
    int rank = H5Sget_simple_extent_dims(space, dims, nullptr);
    hsize_t total = 1;
    for (int i = 0; i < rank; ++i) total *= dims[i];
    out.resize(total);
    herr_t status = H5Dread(ds, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
    H5Sclose(space);
    H5Dclose(ds);
    return status >= 0;
}

static bool read_dataset_i64(hid_t loc_id, const char* name,
                             std::vector<std::int64_t>& out) {
    hid_t ds = H5Dopen2(loc_id, name, H5P_DEFAULT);
    if (ds < 0) return false;
    hid_t space = H5Dget_space(ds);
    hsize_t dims[2];
    int rank = H5Sget_simple_extent_dims(space, dims, nullptr);
    hsize_t total = 1;
    for (int i = 0; i < rank; ++i) total *= dims[i];
    out.resize(total);
    herr_t status = H5Dread(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
    H5Sclose(space);
    H5Dclose(ds);
    return status >= 0;
}

static bool read_attr_str(hid_t loc_id, const char* name, std::string& out) {
    hid_t attr = H5Aopen(loc_id, name, H5P_DEFAULT);
    if (attr < 0) return false;

    hid_t atype = H5Aget_type(attr);
    if (atype < 0) {
        H5Aclose(attr);
        return false;
    }

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

static std::vector<std::string> parse_patch_metadata(const std::string& patches_str) {
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

bool read_mesh_hdf5(const std::string& filename, cfdx::core::Mesh& mesh) {
    mesh.clear();

    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return false;

    // Points (n_points x 3).
    std::vector<double> pts;
    if (!read_dataset_double(file, "points", pts)) {
        pts.clear();
    }
    const std::size_t n_points = pts.size() / 3;
    mesh.points().resize(n_points);
    for (std::size_t i = 0; i < n_points; ++i) {
        mesh.points().set(i, pts[i * 3 + 0], pts[i * 3 + 1], pts[i * 3 + 2]);
    }

    // Face vertices + offsets (64-bit).
    std::vector<std::uint64_t> fv;
    if (!read_dataset_u64(file, "face_vertices", fv)) {
        fv.clear();
    }
    std::vector<std::uint64_t> fo;
    if (!read_dataset_u64(file, "face_offsets", fo)) {
        fo.clear();
    }
    const std::size_t n_faces = fo.size() > 0 ? fo.size() - 1 : 0;
    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::uint64_t off = fo[f];
        const std::uint64_t n = fo[f + 1] - off;
        mesh.faces().push_face(
            std::vector<cfdx::core::FaceConnectivity::Index>(
                fv.begin() + static_cast<std::vector<std::uint64_t>::difference_type>(off), 
                fv.begin() + static_cast<std::vector<std::uint64_t>::difference_type>(off + n)));
    }

    // Owner / neighbour (64-bit).
    std::vector<std::uint64_t> owner;
    if (!read_dataset_u64(file, "owner", owner)) {
        owner.clear();
    }
    std::vector<std::int64_t> neighbour;
    if (!read_dataset_i64(file, "neighbour", neighbour)) {
        neighbour.clear();
    }
    mesh.ownership().resize(owner.size());
    for (std::size_t i = 0; i < owner.size(); ++i) {
        mesh.ownership().set_owner(i, owner[i]);
        mesh.ownership().set_neighbour(i, neighbour[i]);
    }

    // Cell faces + offsets (64-bit).
    std::vector<std::uint64_t> cf;
    if (!read_dataset_u64(file, "cell_faces", cf)) {
        cf.clear();
    }
    std::vector<std::uint64_t> co;
    if (!read_dataset_u64(file, "cell_offsets", co)) {
        co.clear();
    }
    const std::size_t n_cells = co.size() > 0 ? co.size() - 1 : 0;
    for (std::size_t c = 0; c < n_cells; ++c) {
        const std::uint64_t off = co[c];
        const std::uint64_t n = co[c + 1] - off;
        mesh.cells().push_cell(
            std::vector<cfdx::core::CellConnectivity::FaceId>(
                cf.begin() + static_cast<std::vector<std::uint64_t>::difference_type>(off), 
                cf.begin() + static_cast<std::vector<std::uint64_t>::difference_type>(off + n)));
    }

    // Boundary patches.
    std::string patches_str;
    if (read_attr_str(file, "boundary_patches", patches_str)) {
        auto patch_entries = parse_patch_metadata(patches_str);
        
        std::vector<std::uint64_t> patch_face_ids;
        std::vector<std::uint64_t> patch_face_offsets;
        
        if (read_dataset_u64(file, "patch_face_ids", patch_face_ids) &&
            read_dataset_u64(file, "patch_face_offsets", patch_face_offsets)) {
            std::size_t n_patches = patch_entries.size();
            if (n_patches == patch_face_offsets.size() - 1) {
                cfdx::core::BoundaryPatches bp;
                for (std::size_t i = 0; i < n_patches; ++i) {
                    const std::string& entry = patch_entries[i];
                    size_t pos1 = entry.find(':');
                    size_t pos2 = entry.find(':', pos1 + 1);
                    size_t pos3 = entry.find(':', pos2 + 1);
                    
                    if (pos1 != std::string::npos && pos2 != std::string::npos && pos3 != std::string::npos) {
                        std::string name = entry.substr(0, pos1);
                        std::uint64_t start = std::stoull(entry.substr(pos1 + 1, pos2 - pos1 - 1));
                        std::uint64_t count = std::stoull(entry.substr(pos2 + 1, pos3 - pos2 - 1));
                        std::uint32_t type = static_cast<std::uint32_t>(std::stoull(entry.substr(pos3 + 1)));
                        
                        cfdx::core::Patch p;
                        p.name = name;
                        p.type = static_cast<cfdx::core::PatchType>(type);
                        p.face_ids.reserve(count);
                        for (std::uint64_t j = 0; j < count; ++j) {
                            if (start + j < patch_face_ids.size()) {
                                p.face_ids.push_back(patch_face_ids[patch_face_offsets[i] + j]);
                            }
                        }
                        bp.add_patch(p);
                    }
                }
                mesh.set_boundary(bp);
            }
        }
    }

    H5Fclose(file);
    return true;
}

bool read_field_hdf5(const std::string& filename,
                      cfdx::core::Field<double, cfdx::core::Location::CELL>& field) {
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return false;

    hid_t grp = H5Gopen2(file, "fields", H5P_DEFAULT);
    if (grp < 0) { H5Fclose(file); return false; }

    std::vector<double> flat_values;
    if (!read_dataset_double(grp, "values", flat_values)) {
        H5Gclose(grp); H5Fclose(file); return false;
    }

    std::string name, unit, dim_str;
    read_attr_str(grp, "name", name);
    read_attr_str(grp, "unit", unit);
    read_attr_str(grp, "dimension", dim_str);

    std::size_t dim = 1;
    try {
        dim = std::stoul(dim_str);
    } catch (...) { dim = 1; }

    const std::size_t n = flat_values.size() / dim;
    field.set_dimension(dim);
    field.resize(n);

    field.metadata().name = name;
    field.metadata().unit = unit;
    field.metadata().dimension = dim;

    // De-interleave flat array (SoA layout: [x0,y0,z0, x1,y1,z1, ...])
    for (std::size_t c = 0; c < dim; ++c) {
        double* comp_data = field.component_data(c);
        for (std::size_t i = 0; i < n; ++i) {
            comp_data[i] = flat_values[i * dim + c];
        }
    }

    H5Gclose(grp);
    H5Fclose(file);
    return true;
}

}  // namespace io
}  // namespace cfdx
