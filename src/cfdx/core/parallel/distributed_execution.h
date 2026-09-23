#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/parallel/mpi_utils.h"
#include "cfdx/core/parallel/mesh_partitioner.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <utility>

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

        // A cell-centred halo is bidirectional across an MPI interface:
        // each rank publishes its local cell value to the other rank.
        if (owner_rank == rank) {
            halo.send_local_cells[static_cast<std::size_t>(neighbour_rank)]
                .push_back(field.local_index(static_cast<std::uint64_t>(owner)));
            halo.recv_global_ids[static_cast<std::size_t>(neighbour_rank)]
                .push_back(static_cast<std::uint64_t>(neighbour));
        }
        if (neighbour_rank == rank) {
            halo.send_local_cells[static_cast<std::size_t>(owner_rank)]
                .push_back(field.local_index(static_cast<std::uint64_t>(neighbour)));
            halo.recv_global_ids[static_cast<std::size_t>(owner_rank)]
                .push_back(static_cast<std::uint64_t>(owner));
        }
    }

    for (auto& v : halo.send_local_cells) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }
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

    // Result ordering is deterministic: peer order, then the halo plan's
    // receive-ID order, with all field components contiguous per cell.
    std::unordered_map<std::uint64_t, std::vector<double>> received;
    const std::size_t dimension_size = field.dimension();
    if (dimension_size == 0)
        throw std::invalid_argument("exchange_distributed_cell_halo: field dimension must be positive");
    if (dimension_size > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error("exchange_distributed_cell_halo: dimension overflow");
    const int dimension = static_cast<int>(dimension_size);

    for (int peer = 0; peer < size; ++peer) {
        if (peer == rank) continue;
        const auto& send = halo.send_local_cells[static_cast<std::size_t>(peer)];
        const auto& recv = halo.recv_global_ids[static_cast<std::size_t>(peer)];
        if (send.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
            recv.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::overflow_error("exchange_distributed_cell_halo: cell count overflow");
        const int send_count = static_cast<int>(send.size());
        const int recv_count = static_cast<int>(recv.size());
        int remote_send_count = 0;
        int remote_recv_count = 0;
        int remote_dimension = 0;
        MPI_Status status{};
        MPI_Sendrecv(&send_count, 1, MPI_INT, peer, 610,
                     &remote_send_count, 1, MPI_INT, peer, 610,
                     comm, &status);
        MPI_Sendrecv(&recv_count, 1, MPI_INT, peer, 611,
                     &remote_recv_count, 1, MPI_INT, peer, 611,
                     comm, &status);
        MPI_Sendrecv(&dimension, 1, MPI_INT, peer, 613,
                     &remote_dimension, 1, MPI_INT, peer, 613,
                     comm, &status);
        if (remote_send_count != recv_count || remote_recv_count != send_count)
            throw std::runtime_error("exchange_distributed_cell_halo: asymmetric peer counts");
        if (remote_dimension != dimension)
            throw std::invalid_argument("exchange_distributed_cell_halo: peer field dimensions differ");

        const std::size_t send_values =
            static_cast<std::size_t>(send_count) * dimension_size;
        const std::size_t recv_values =
            static_cast<std::size_t>(recv_count) * dimension_size;
        if (send_values > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
            recv_values > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::overflow_error("exchange_distributed_cell_halo: MPI value count overflow");

        std::vector<double> send_buffer(send_values);
        for (int i = 0; i < send_count; ++i)
            for (std::size_t comp = 0; comp < dimension_size; ++comp)
                send_buffer[static_cast<std::size_t>(i) * dimension_size + comp] =
                    field(send[static_cast<std::size_t>(i)], comp);

        std::vector<double> recv_buffer(recv_values);
        MPI_Sendrecv(send_buffer.empty() ? nullptr : send_buffer.data(),
                     static_cast<int>(send_values), MPI_DOUBLE, peer, 612,
                     recv_buffer.empty() ? nullptr : recv_buffer.data(),
                     static_cast<int>(recv_values), MPI_DOUBLE, peer, 612,
                     comm, &status);

        for (int i = 0; i < recv_count; ++i) {
            std::vector<double> values(dimension_size);
            for (std::size_t comp = 0; comp < dimension_size; ++comp)
                values[comp] =
                    recv_buffer[static_cast<std::size_t>(i) * dimension_size + comp];
            const auto [it, inserted] =
                received.emplace(recv[static_cast<std::size_t>(i)], std::move(values));
            if (!inserted)
                throw std::runtime_error("exchange_distributed_cell_halo: duplicate received global id");
        }
    }

    std::vector<double> result;
    for (const auto& ids_for_peer : halo.recv_global_ids) {
        for (const auto gid : ids_for_peer) {
            const auto it = received.find(gid);
            if (it == received.end())
                throw std::runtime_error("exchange_distributed_cell_halo: missing received global id");
            result.insert(result.end(), it->second.begin(), it->second.end());
        }
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

    // Exact distributed uniqueness check without replicating the complete
    // ownership list on every rank. IDs are routed to an owner rank, where
    // duplicates are detected locally.
    const int size = mpi_size(comm);
    std::vector<std::vector<std::uint64_t>> buckets(static_cast<std::size_t>(size));
    for (const auto id : local_ids)
        buckets[static_cast<std::size_t>(id % static_cast<std::uint64_t>(size))].push_back(id);

    std::vector<int> send_counts(static_cast<std::size_t>(size), 0);
    std::vector<int> recv_counts(static_cast<std::size_t>(size), 0);
    std::vector<int> send_displs(static_cast<std::size_t>(size), 0);
    std::vector<int> recv_displs(static_cast<std::size_t>(size), 0);
    for (int r = 0; r < size; ++r) {
        if (buckets[static_cast<std::size_t>(r)].size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::overflow_error("validate_distributed_ids: MPI count overflow");
        send_counts[static_cast<std::size_t>(r)] =
            static_cast<int>(buckets[static_cast<std::size_t>(r)].size());
        if (r > 0)
            send_displs[static_cast<std::size_t>(r)] =
                send_displs[static_cast<std::size_t>(r - 1)] +
                send_counts[static_cast<std::size_t>(r - 1)];
    }
    std::vector<std::uint64_t> send_buffer(local_ids.size());
    for (int r = 0; r < size; ++r)
        std::copy(buckets[static_cast<std::size_t>(r)].begin(),
                  buckets[static_cast<std::size_t>(r)].end(),
                  send_buffer.begin() + send_displs[static_cast<std::size_t>(r)]);

    MPI_Alltoall(send_counts.data(), 1, MPI_INT,
                 recv_counts.data(), 1, MPI_INT, comm);
    for (int r = 1; r < size; ++r)
        recv_displs[static_cast<std::size_t>(r)] =
            recv_displs[static_cast<std::size_t>(r - 1)] +
            recv_counts[static_cast<std::size_t>(r - 1)];
    const std::size_t recv_total = static_cast<std::size_t>(
        recv_displs.back() + recv_counts.back());
    std::vector<std::uint64_t> received(recv_total);
    MPI_Alltoallv(send_buffer.empty() ? nullptr : send_buffer.data(),
                  send_counts.data(), send_displs.data(), MPI_UINT64_T,
                  received.empty() ? nullptr : received.data(),
                  recv_counts.data(), recv_displs.data(), MPI_UINT64_T, comm);

    std::sort(received.begin(), received.end());
    if (std::adjacent_find(received.begin(), received.end()) != received.end())
        throw std::invalid_argument("validate_distributed_ids: duplicate global id");

    const std::int64_t unique_count =
        static_cast<std::int64_t>(received.size());
    const std::int64_t global_unique = mpi_allreduce_sum(unique_count, comm);
    if (global_unique != static_cast<std::int64_t>(global_size))
        throw std::invalid_argument(
            "validate_distributed_ids: ownership does not cover global domain");
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
    const hsize_t dims[2] = {
        static_cast<hsize_t>(state.dimension()),
        static_cast<hsize_t>(state.global_size())
    };
    hid_t fs = H5Screate_simple(2, dims, nullptr);
    if (fs < 0) throw std::runtime_error("distributed HDF5: failed to create field dataspace");
    hid_t ds = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, fs,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(fs);
    if (ds < 0) throw std::runtime_error(std::string("distributed HDF5: failed to create ") + name);

    fs = H5Dget_space(ds);
    const hsize_t mem_dims[2] = {
        static_cast<hsize_t>(state.dimension()),
        static_cast<hsize_t>(state.local_size())
    };
    hid_t ms = H5Screate_simple(2, mem_dims, nullptr);
    if (ms < 0) { H5Sclose(fs); H5Dclose(ds); throw std::runtime_error("distributed HDF5: failed to create field memory dataspace"); }

    std::vector<hsize_t> coords(2 * state.local_size() * state.dimension());
    for (std::size_t comp = 0; comp < state.dimension(); ++comp) {
        for (std::size_t i = 0; i < state.local_size(); ++i) {
            const std::size_t k = comp * state.local_size() + i;
            coords[2 * k] = static_cast<hsize_t>(comp);
            coords[2 * k + 1] = static_cast<hsize_t>(state.global_ids()[i]);
        }
    }
    if (coords.empty()) {
        if (H5Sselect_none(fs) < 0) {
            H5Sclose(ms); H5Sclose(fs); H5Dclose(ds);
            throw std::runtime_error("distributed HDF5: failed to select empty field");
        }
    } else if (H5Sselect_elements(fs, H5S_SELECT_SET, state.local_size() * state.dimension(),
                                  coords.data()) < 0) {
        H5Sclose(ms); H5Sclose(fs); H5Dclose(ds);
        throw std::runtime_error("distributed HDF5: failed to select global cell ids");
    }

    std::vector<double> values(state.local_size() * state.dimension());
    for (std::size_t comp = 0; comp < state.dimension(); ++comp)
        for (std::size_t i = 0; i < state.local_size(); ++i)
            values[comp * state.local_size() + i] = state(i, comp);

    hid_t dxpl = collective_dxpl();
    const herr_t status = H5Dwrite(ds, H5T_NATIVE_DOUBLE, ms, fs, dxpl,
                                   values.empty() ? nullptr : values.data());
    H5Pclose(dxpl);
    H5Sclose(ms);
    H5Sclose(fs);
    H5Dclose(ds);
    if (status < 0)
        throw std::runtime_error(std::string("distributed HDF5: failed to write ") + name);
}

inline void read_dataset(hid_t file, const char* name,
                         DistributedCellField& state)
{
    hid_t ds = H5Dopen2(file, name, H5P_DEFAULT);
    if (ds < 0) throw std::runtime_error(std::string("distributed HDF5: missing ") + name);
    hid_t fs = H5Dget_space(ds);
    int rank = H5Sget_simple_extent_ndims(fs);
    hsize_t dims[2] = {0, 0};
    if (rank != 2 || H5Sget_simple_extent_dims(fs, dims, nullptr) < 0 ||
        dims[0] != static_cast<hsize_t>(state.dimension()) ||
        dims[1] != static_cast<hsize_t>(state.global_size())) {
        H5Sclose(fs); H5Dclose(ds);
        throw std::runtime_error("distributed HDF5: field dataset shape mismatch");
    }

    const hsize_t mem_dims[2] = {
        static_cast<hsize_t>(state.dimension()),
        static_cast<hsize_t>(state.local_size())
    };
    hid_t ms = H5Screate_simple(2, mem_dims, nullptr);
    if (ms < 0) { H5Sclose(fs); H5Dclose(ds); throw std::runtime_error("distributed HDF5: failed to create field memory dataspace"); }

    std::vector<hsize_t> coords(2 * state.local_size() * state.dimension());
    for (std::size_t comp = 0; comp < state.dimension(); ++comp)
        for (std::size_t i = 0; i < state.local_size(); ++i) {
            const std::size_t k = comp * state.local_size() + i;
            coords[2 * k] = static_cast<hsize_t>(comp);
            coords[2 * k + 1] = static_cast<hsize_t>(state.global_ids()[i]);
        }
    if (coords.empty()) {
        if (H5Sselect_none(fs) < 0) {
            H5Sclose(ms); H5Sclose(fs); H5Dclose(ds);
            throw std::runtime_error("distributed HDF5: failed to select empty field");
        }
    } else if (H5Sselect_elements(fs, H5S_SELECT_SET, state.local_size() * state.dimension(),
                                  coords.data()) < 0) {
        H5Sclose(ms); H5Sclose(fs); H5Dclose(ds);
        throw std::runtime_error("distributed HDF5: failed to select global cell ids");
    }

    std::vector<double> values(state.local_size() * state.dimension());
    hid_t dxpl = collective_dxpl();
    const herr_t status = H5Dread(ds, H5T_NATIVE_DOUBLE, ms, fs, dxpl,
                                  values.empty() ? nullptr : values.data());
    H5Pclose(dxpl);
    H5Sclose(ms);
    H5Sclose(fs);
    H5Dclose(ds);
    if (status < 0)
        throw std::runtime_error(std::string("distributed HDF5: failed to read ") + name);
    for (std::size_t comp = 0; comp < state.dimension(); ++comp)
        for (std::size_t i = 0; i < state.local_size(); ++i)
            state(i, comp) = values[comp * state.local_size() + i];
}

} // namespace detail

inline void write_distributed_checkpoint(
    const std::string& path,
    const DistributedCellField& state,
    MPI_Comm comm = MPI_COMM_WORLD)
{
    validate_distributed_ids(state.global_ids(), state.global_size(), comm);
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

    detail::write_dataset(file, "field", state);
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
    // Validate the checkpoint identity at the target rank's selected IDs.
    // This makes the persistent global-cell-ID dataset part of the restart
    // contract rather than merely redundant metadata.
    hid_t ids_verify_fs = H5Dget_space(ids_ds);
    hid_t ids_verify_ms = detail::local_memspace(state.local_size());
    detail::select_ids(ids_verify_fs, state.global_ids());
    hid_t ids_verify_dxpl = detail::collective_dxpl();
    std::vector<std::uint64_t> checkpoint_ids(state.local_size());
    if (H5Dread(ids_ds, H5T_NATIVE_UINT64, ids_verify_ms, ids_verify_fs,
                ids_verify_dxpl, checkpoint_ids.data()) < 0) {
        H5Pclose(ids_verify_dxpl); H5Sclose(ids_verify_ms);
        H5Sclose(ids_verify_fs); H5Sclose(ids_fs); H5Dclose(ids_ds); H5Fclose(file);
        throw std::runtime_error("distributed HDF5: failed to read global_cell_ids");
    }
    H5Pclose(ids_verify_dxpl);
    H5Sclose(ids_verify_ms);
    H5Sclose(ids_verify_fs);
    for (std::size_t i = 0; i < checkpoint_ids.size(); ++i) {
        if (checkpoint_ids[i] != state.global_ids()[i]) {
            H5Sclose(ids_fs); H5Dclose(ids_ds); H5Fclose(file);
            throw std::runtime_error("distributed HDF5: checkpoint global cell ID mismatch");
        }
    }
    H5Sclose(ids_fs);
    H5Dclose(ids_ds);

    detail::read_dataset(file, "field", state);
    H5Fclose(file);
}

#endif // CFDX_ENABLE_PARALLEL_HDF5

} // namespace cfdx::core::parallel
