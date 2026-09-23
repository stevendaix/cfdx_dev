#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/parallel/mpi_utils.h"
#include "cfdx/core/parallel/mesh_partitioner.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#ifdef CFDX_ENABLE_PARALLEL_HDF5
#include <H5Dpublic.h>
#include <H5Fpublic.h>
#include <H5Ppublic.h>
#include <H5Spublic.h>
#include <H5Tpublic.h>
#include <H5FDmpi.h>
#include <H5FDmpio.h>
#endif

namespace cfdx::core::parallel {

/*
 * Phase 5 distributed execution state.
 *
 * A DistributedField owns only the cells assigned to this MPI rank.  The
 * persistent global cell id is the identity used by halo exchange and
 * checkpoints; local vector order is deliberately not part of the restart
 * contract.
 */
class DistributedCellField {
public:
    DistributedCellField() = default;

    DistributedCellField(std::size_t global_size,
                         std::vector<std::uint64_t> global_ids,
                         std::size_t dimension = 1,
                         const std::string& name = "")
        : global_size_(global_size),
          global_ids_(std::move(global_ids)),
          field_(global_ids_.size(), name, "", dimension)
    {
        validate_ids();
    }

    std::size_t global_size() const noexcept { return global_size_; }
    std::size_t local_size() const noexcept { return global_ids_.size(); }
    std::size_t dimension() const noexcept { return field_.dimension(); }

    const std::vector<std::uint64_t>& global_ids() const noexcept { return global_ids_; }

    std::size_t local_index(std::uint64_t global_id) const {
        const auto it = local_index_.find(global_id);
        if (it == local_index_.end())
            throw std::out_of_range("DistributedCellField: global cell id is not local");
        return it->second;
    }

    cfdx::core::Field<double, cfdx::core::Location::CELL>& field() noexcept {
        return field_;
    }
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& field() const noexcept {
        return field_;
    }

    double& operator()(std::size_t local, std::size_t component = 0) {
        return field_(local, component);
    }
    double operator()(std::size_t local, std::size_t component = 0) const {
        return field_(local, component);
    }

private:
    void validate_ids() {
        if (global_size_ == 0 && !global_ids_.empty())
            throw std::invalid_argument("DistributedCellField: non-empty ids with zero global size");
        std::unordered_set<std::uint64_t> seen;
        seen.reserve(global_ids_.size());
        local_index_.reserve(global_ids_.size());
        for (std::size_t i = 0; i < global_ids_.size(); ++i) {
            if (global_ids_[i] >= global_size_)
                throw std::invalid_argument("DistributedCellField: global id out of range");
            if (!seen.insert(global_ids_[i]).second)
                throw std::invalid_argument("DistributedCellField: duplicate local global id");
            local_index_.emplace(global_ids_[i], i);
        }
    }

    std::size_t global_size_ = 0;
    std::vector<std::uint64_t> global_ids_;
    std::unordered_map<std::uint64_t, std::size_t> local_index_;
    cfdx::core::Field<double, cfdx::core::Location::CELL> field_;
};

inline std::vector<std::uint64_t> owned_global_cell_ids(
    const Partition& partition, int rank)
{
    if (rank < 0 || rank >= partition.n_parts)
        throw std::invalid_argument("owned_global_cell_ids: invalid rank");
    std::vector<std::uint64_t> ids;
    for (std::size_t c = 0; c < partition.cell_rank.size(); ++c) {
        if (partition.cell_rank[c] == rank)
            ids.push_back(static_cast<std::uint64_t>(c));
    }
    return ids;
}

struct DistributedHalo {
    std::vector<std::vector<std::size_t>> send_local_cells;
    std::vector<std::vector<std::uint64_t>> recv_global_ids;
};

inline DistributedHalo build_distributed_halo(
    const Mesh& mesh,
    const Partition& partition,
    const DistributedCellField& field,
    MPI_Comm comm = MPI_COMM_WORLD)
{
    const int rank = mpi_rank(comm);
    const int size = mpi_size(comm);
    if (partition.n_parts != size)
        throw std::invalid_argument("build_distributed_halo: partition size differs from communicator");

    DistributedHalo halo;
    halo.send_local_cells.resize(static_cast<std::size_t>(size));
    halo.recv_global_ids.resize(static_cast<std::size_t>(size));

    const auto& own = mesh.ownership();
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const auto owner = own.owner(f);
        const auto neighbour = own.neighbour(f);
        if (neighbour < 0) continue;
        const int owner_rank = partition.cell_rank[owner];
        const int neighbour_rank = partition.cell_rank[static_cast<std::size_t>(neighbour)];
        if (owner_rank == neighbour_rank) continue;

        if (owner_rank == rank) {
            halo.send_local_cells[static_cast<std::size_t>(neighbour_rank)]
                .push_back(field.local_index(static_cast<std::uint64_t>(owner)));
        }
        if (neighbour_rank == rank) {
            halo.recv_global_ids[static_cast<std::size_t>(owner_rank)]
                .push_back(static_cast<std::uint64_t>(owner));
        }
    }

    for (auto& v : halo.send_local_cells)
        std::sort(v.begin(), v.end());
    for (auto& v : halo.recv_global_ids) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
    return halo;
}

inline std::vector<double> exchange_distributed_cell_halo(
    const DistributedCellField& field,
    const DistributedHalo& halo,
    MPI_Comm comm = MPI_COMM_WORLD)
{
    const int rank = mpi_rank(comm);
    const int size = mpi_size(comm);
    if (static_cast<int>(halo.send_local_cells.size()) != size ||
        static_cast<int>(halo.recv_global_ids.size()) != size)
        throw std::invalid_argument("exchange_distributed_cell_halo: invalid halo plan");

    std::unordered_map<std::uint64_t, double> received;
    for (int peer = 0; peer < size; ++peer) {
        if (peer == rank) continue;
        const auto& send = halo.send_local_cells[static_cast<std::size_t>(peer)];
        const auto& recv = halo.recv_global_ids[static_cast<std::size_t>(peer)];
        const int send_count = static_cast<int>(send.size());
        const int recv_count = static_cast<int>(recv.size());
        int remote_send_count = 0;
        int remote_recv_count = 0;
        MPI_Status status{};
        MPI_Sendrecv(&send_count, 1, MPI_INT, peer, 610,
                     &remote_send_count, 1, MPI_INT, peer, 610,
                     comm, &status);
        MPI_Sendrecv(&recv_count, 1, MPI_INT, peer, 611,
                     &remote_recv_count, 1, MPI_INT, peer, 611,
                     comm, &status);
        if (remote_send_count != recv_count || remote_recv_count != send_count)
            throw std::runtime_error("exchange_distributed_cell_halo: asymmetric peer counts");

        std::vector<double> send_buffer(send_count);
        for (int i = 0; i < send_count; ++i)
            send_buffer[static_cast<std::size_t>(i)] = field(send[static_cast<std::size_t>(i)]);
        std::vector<double> recv_buffer(recv_count);
        MPI_Sendrecv(send_buffer.empty() ? nullptr : send_buffer.data(), send_count,
                     MPI_DOUBLE, peer, 612,
                     recv_buffer.empty() ? nullptr : recv_buffer.data(), recv_count,
                     MPI_DOUBLE, peer, 612, comm, &status);
        for (int i = 0; i < recv_count; ++i)
            received.emplace(recv[static_cast<std::size_t>(i)],
                             recv_buffer[static_cast<std::size_t>(i)]);
    }

    std::vector<double> result;
    result.reserve(received.size());
    for (const auto& [id, value] : received) {
        (void)id;
        result.push_back(value);
    }
    return result;
}

inline void validate_distributed_ids(
    const std::vector<std::uint64_t>& local_ids,
    std::size_t global_size,
    MPI_Comm comm = MPI_COMM_WORLD)
{
    std::unordered_set<std::uint64_t> local;
    local.reserve(local_ids.size());
    for (const auto id : local_ids) {
        if (id >= global_size)
            throw std::invalid_argument("validate_distributed_ids: id out of range");
        if (!local.insert(id).second)
            throw std::invalid_argument("validate_distributed_ids: duplicate local id");
    }

    const auto gathered = mpi_allgather(local_ids, comm);
    if (gathered.size() != global_size)
        throw std::invalid_argument("validate_distributed_ids: global ownership count mismatch");

    std::unordered_set<std::uint64_t> global;
    global.reserve(gathered.size());
    for (const auto id : gathered) {
        if (!global.insert(id).second)
            throw std::invalid_argument("validate_distributed_ids: duplicate global id");
    }
    if (global.size() != global_size)
        throw std::invalid_argument("validate_distributed_ids: ownership does not cover global domain");
}

#ifdef CFDX_ENABLE_PARALLEL_HDF5

namespace detail {

inline hid_t distributed_fapl(MPI_Comm comm) {
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl < 0 || H5Pset_fapl_mpio(fapl, comm, MPI_INFO_NULL) < 0) {
        if (fapl >= 0) H5Pclose(fapl);
        throw std::runtime_error("distributed HDF5: failed to create MPI file access property list");
    }
    return fapl;
}

inline hid_t collective_dxpl() {
    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
    if (dxpl < 0 || H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE) < 0) {
        if (dxpl >= 0) H5Pclose(dxpl);
        throw std::runtime_error("distributed HDF5: failed to create collective transfer property list");
    }
    return dxpl;
}

inline void select_ids(hid_t dataspace,
                       const std::vector<std::uint64_t>& ids)
{
    if (ids.empty()) {
        if (H5Sselect_none(dataspace) < 0)
            throw std::runtime_error("distributed HDF5: H5Sselect_none failed");
        return;
    }
    std::vector<hsize_t> coords(ids.size());
    for (std::size_t i = 0; i < ids.size(); ++i)
        coords[i] = static_cast<hsize_t>(ids[i]);
    if (H5Sselect_elements(dataspace, H5S_SELECT_SET, coords.size(), coords.data()) < 0)
        throw std::runtime_error("distributed HDF5: failed to select global cell ids");
}

inline hid_t local_memspace(std::size_t n) {
    const hsize_t count = static_cast<hsize_t>(n);
    hid_t space = H5Screate_simple(1, &count, nullptr);
    if (space < 0) throw std::runtime_error("distributed HDF5: failed to create memory dataspace");
    return space;
}

inline void write_dataset(hid_t file, const char* name,
                          const DistributedCellField& state)
{
    const hsize_t global_n = static_cast<hsize_t>(state.global_size());
    hid_t fs = H5Screate_simple(1, &global_n, nullptr);
    if (fs < 0) throw std::runtime_error("distributed HDF5: failed to create file dataspace");
    hid_t ds = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, fs,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(fs);
    if (ds < 0) throw std::runtime_error(std::string("distributed HDF5: failed to create ") + name);

    fs = H5Dget_space(ds);
    hid_t ms = local_memspace(state.local_size());
    select_ids(fs, state.global_ids());
    hid_t dxpl = collective_dxpl();
    for (std::size_t c = 0; c < state.dimension(); ++c) {
        std::vector<double> values(state.local_size());
        for (std::size_t i = 0; i < state.local_size(); ++i)
            values[i] = state(i, c);
        if (H5Dwrite(ds, H5T_NATIVE_DOUBLE, ms, fs, dxpl, values.data()) < 0) {
            H5Pclose(dxpl); H5Sclose(ms); H5Sclose(fs); H5Dclose(ds);
            throw std::runtime_error(std::string("distributed HDF5: failed to write ") + name);
        }
    }
    H5Pclose(dxpl);
    H5Sclose(ms);
    H5Sclose(fs);
    H5Dclose(ds);
}

inline void read_dataset(hid_t file, const char* name,
                         DistributedCellField& state)
{
    hid_t ds = H5Dopen2(file, name, H5P_DEFAULT);
    if (ds < 0) throw std::runtime_error(std::string("distributed HDF5: missing ") + name);
    hid_t fs = H5Dget_space(ds);
    int rank = H5Sget_simple_extent_ndims(fs);
    if (rank != 1) {
        H5Sclose(fs); H5Dclose(ds);
        throw std::runtime_error("distributed HDF5: field dataset must be one-dimensional");
    }
    hsize_t dims[1] = {0};
    H5Sget_simple_extent_dims(fs, dims, nullptr);
    if (dims[0] != static_cast<hsize_t>(state.global_size())) {
        H5Sclose(fs); H5Dclose(ds);
        throw std::runtime_error("distributed HDF5: global size mismatch");
    }

    hid_t ms = local_memspace(state.local_size());
    select_ids(fs, state.global_ids());
    hid_t dxpl = collective_dxpl();
    std::vector<double> values(state.local_size());
    for (std::size_t c = 0; c < state.dimension(); ++c) {
        if (H5Dread(ds, H5T_NATIVE_DOUBLE, ms, fs, dxpl, values.data()) < 0) {
            H5Pclose(dxpl); H5Sclose(ms); H5Sclose(fs); H5Dclose(ds);
            throw std::runtime_error(std::string("distributed HDF5: failed to read ") + name);
        }
        for (std::size_t i = 0; i < state.local_size(); ++i)
            state(i, c) = values[i];
    }
    H5Pclose(dxpl);
    H5Sclose(ms);
    H5Sclose(fs);
    H5Dclose(ds);
}

} // namespace detail

inline void write_distributed_checkpoint(
    const std::string& path,
    const DistributedCellField& state,
    MPI_Comm comm = MPI_COMM_WORLD)
{
    validate_distributed_ids(state.global_ids(), state.global_size(), comm);
    const int rank = mpi_rank(comm);
    hid_t fapl = detail::distributed_fapl(comm);
    hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    H5Pclose(fapl);
    if (file < 0) throw std::runtime_error("distributed HDF5: failed to create checkpoint");

    const hsize_t global_n = static_cast<hsize_t>(state.global_size());
    hid_t ids_space = H5Screate_simple(1, &global_n, nullptr);
    hid_t ids_ds = H5Dcreate2(file, "global_cell_ids", H5T_NATIVE_UINT64,
                              ids_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(ids_space);
    if (ids_ds < 0) { H5Fclose(file); throw std::runtime_error("distributed HDF5: failed to create global_cell_ids"); }

    hid_t ids_fs = H5Dget_space(ids_ds);
    hid_t ids_ms = detail::local_memspace(state.local_size());
    detail::select_ids(ids_fs, state.global_ids());
    hid_t ids_dxpl = detail::collective_dxpl();
    std::vector<std::uint64_t> ids = state.global_ids();
    if (H5Dwrite(ids_ds, H5T_NATIVE_UINT64, ids_ms, ids_fs, ids_dxpl, ids.data()) < 0) {
        H5Pclose(ids_dxpl); H5Sclose(ids_ms); H5Sclose(ids_fs); H5Dclose(ids_ds); H5Fclose(file);
        throw std::runtime_error("distributed HDF5: failed to write global_cell_ids");
    }
    H5Pclose(ids_dxpl); H5Sclose(ids_ms); H5Sclose(ids_fs); H5Dclose(ids_ds);

    write_dataset(file, "field", state);
    if (rank == 0) {
        // A scalar marker makes the on-disk contract self-describing.
        hid_t attr_space = H5Screate(H5S_SCALAR);
        hid_t attr = H5Acreate2(file, "cfdx_distributed_checkpoint_version",
                                H5T_NATIVE_INT, attr_space, H5P_DEFAULT, H5P_DEFAULT);
        const int version = 1;
        H5Awrite(attr, H5T_NATIVE_INT, &version);
        H5Aclose(attr);
        H5Sclose(attr_space);
    }
    MPI_Barrier(comm);
    H5Fclose(file);
}

inline void read_distributed_checkpoint(
    const std::string& path,
    DistributedCellField& state,
    MPI_Comm comm = MPI_COMM_WORLD)
{
    validate_distributed_ids(state.global_ids(), state.global_size(), comm);
    hid_t fapl = detail::distributed_fapl(comm);
    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, fapl);
    H5Pclose(fapl);
    if (file < 0) throw std::runtime_error("distributed HDF5: failed to open checkpoint");

    hid_t ids_ds = H5Dopen2(file, "global_cell_ids", H5P_DEFAULT);
    if (ids_ds < 0) { H5Fclose(file); throw std::runtime_error("distributed HDF5: missing global_cell_ids"); }
    hid_t ids_fs = H5Dget_space(ids_ds);
    hsize_t dims[1] = {0};
    H5Sget_simple_extent_dims(ids_fs, dims, nullptr);
    if (dims[0] != static_cast<hsize_t>(state.global_size())) {
        H5Sclose(ids_fs); H5Dclose(ids_ds); H5Fclose(file);
        throw std::runtime_error("distributed HDF5: checkpoint global size mismatch");
    }
    H5Sclose(ids_fs);
    H5Dclose(ids_ds);

    detail::read_dataset(file, "field", state);
    H5Fclose(file);
}

#endif // CFDX_ENABLE_PARALLEL_HDF5

} // namespace cfdx::core::parallel
