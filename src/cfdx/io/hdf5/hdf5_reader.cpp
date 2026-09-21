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
#include <sstream>
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
    herr_t status = total == 0 ? 0 : H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
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
    herr_t status = total == 0 ? 0 : H5Dread(ds, H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
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
    herr_t status = total == 0 ? 0 : H5Dread(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
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
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return false;

    auto fail = [&](const std::string& reason) {
        std::cerr << "HDF5 integrity error in '" << filename
                  << "': " << reason << '\n';
        H5Fclose(file);
        return false;
    };

    std::vector<double> pts;
    std::vector<std::uint64_t> fv, fo, owner, cf, co;
    std::vector<std::int64_t> neighbour;

    // A CFDX mesh is valid only if all core topology datasets are present.
    if (!read_dataset_double(file, "points", pts) ||
        !read_dataset_u64(file, "face_vertices", fv) ||
        !read_dataset_u64(file, "face_offsets", fo) ||
        !read_dataset_u64(file, "owner", owner) ||
        !read_dataset_i64(file, "neighbour", neighbour) ||
        !read_dataset_u64(file, "cell_faces", cf) ||
        !read_dataset_u64(file, "cell_offsets", co)) {
        return fail("missing required topology dataset");
    }

    if (pts.size() % 3 != 0 || fo.empty() || co.empty() ||
        fo.front() != 0 || co.front() != 0 ||
        fo.back() != fv.size() || co.back() != cf.size() ||
        owner.size() != neighbour.size() ||
        fo.size() - 1 != owner.size()) {
        return fail("invalid topology dataset dimensions or CSR terminal offsets");
    }

    const std::size_t n_points = pts.size() / 3;
    const std::size_t n_faces = fo.size() - 1;
    const std::size_t n_cells = co.size() - 1;

    // Validate CSR offsets before converting uint64_t to size_t.
    for (std::size_t i = 1; i < fo.size(); ++i) {
        if (fo[i] < fo[i - 1] || fo[i] > fv.size()) return fail("face CSR offsets are not monotonic or exceed face-vertex storage");
    }
    for (std::size_t i = 1; i < co.size(); ++i) {
        if (co[i] < co[i - 1] || co[i] > cf.size()) return fail("cell CSR offsets are not monotonic or exceed cell-face storage");
    }

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
            if (fv[j] >= n_points) return fail("face vertex index is outside the point array");
            vertices.push_back(
                static_cast<cfdx::core::FaceConnectivity::Index>(fv[j]));
        }
        if (vertices.size() < 3) return fail("face contains fewer than three vertices");
        mesh.faces().push_face(std::move(vertices));
    }

    mesh.ownership().resize(n_faces);
    for (std::size_t i = 0; i < n_faces; ++i) {
        if (owner[i] >= n_cells) return fail("owner index is outside the cell range");
        if (neighbour[i] < -1 ||
            (neighbour[i] >= 0 &&
             static_cast<std::uint64_t>(neighbour[i]) >= n_cells)) {
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
            if (cf[j] >= n_faces) return fail("cell-face index is outside the face range");
            faces.push_back(
                static_cast<cfdx::core::CellConnectivity::FaceId>(cf[j]));
        }
        if (faces.empty()) return fail("cell contains no faces");
        mesh.cells().push_cell(std::move(faces));
    }

    // Boundary patches are optional for an otherwise valid topological mesh,
    // but when present their metadata and CSR representation must agree.
    std::string patches_str;
    if (read_attr_str(file, "boundary_patches", patches_str)) {
        const auto patch_entries = parse_patch_metadata(patches_str);
        std::vector<std::uint64_t> patch_face_ids, patch_face_offsets;
        if (!read_dataset_u64(file, "patch_face_ids", patch_face_ids) ||
            !read_dataset_u64(file, "patch_face_offsets", patch_face_offsets) ||
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
            if (p1 == std::string::npos || p2 == std::string::npos ||
                p3 == std::string::npos) {
                return fail("boundary patch metadata entry is malformed");
            }

            cfdx::core::Patch patch;
            patch.name = entry.substr(0, p1);
            try {
                const std::uint64_t start = std::stoull(
                    entry.substr(p1 + 1, p2 - p1 - 1));
                const std::uint64_t count = std::stoull(
                    entry.substr(p2 + 1, p3 - p2 - 1));
                patch.type = static_cast<cfdx::core::PatchType>(
                    std::stoull(entry.substr(p3 + 1)));

                const std::uint64_t offset_begin = patch_face_offsets[i];
                const std::uint64_t offset_end = patch_face_offsets[i + 1];
                if (count != offset_end - offset_begin ||
                    start > n_faces ||
                    start + count > n_faces) {
                    return fail("boundary patch face count/range is inconsistent");
                }

                patch.face_ids.reserve(static_cast<std::size_t>(count));
                for (std::uint64_t j = 0; j < count; ++j) {
                    const std::uint64_t face_id = patch_face_ids[offset_begin + j];
                    if (face_id >= n_faces) return fail("boundary patch references a face outside the mesh");
                    patch.face_ids.push_back(
                        static_cast<cfdx::core::FaceIndex>(face_id));
                }
            } catch (...) {
                return fail("mesh topology validation failed");
            }
            bp.add_patch(std::move(patch));
        }
        mesh.set_boundary(bp);
    }

    const bool valid = mesh.topo_validate().ok;
    H5Fclose(file);
    return valid;
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