#include <H5public.h>
#include <H5Dpublic.h>
#include <H5Fpublic.h>
#include <H5Ppublic.h>
#include <H5FDmpi.h>
#include <H5FDmpio.h>
#include <H5Spublic.h>
#include <mpi.h>

#include <cassert>
#include <cstdio>
#include <cstdint>
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    assert(size == 2);

    constexpr hsize_t local_n = 4;
    const hsize_t global_n = local_n * static_cast<hsize_t>(size);
    const char* filename = "parallel_hdf5_hyperslab_test.h5";

    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    assert(fapl >= 0);
    assert(H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL) >= 0);

    hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    assert(file >= 0);
    H5Pclose(fapl);

    hid_t filespace = H5Screate_simple(1, &global_n, nullptr);
    assert(filespace >= 0);
    hid_t dataset = H5Dcreate2(file, "values", H5T_NATIVE_DOUBLE, filespace,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(dataset >= 0);

    const hsize_t offset = static_cast<hsize_t>(rank) * local_n;
    const hsize_t count = local_n;
    assert(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, nullptr,
                               &count, nullptr) >= 0);

    hid_t memspace = H5Screate_simple(1, &count, nullptr);
    assert(memspace >= 0);
    std::vector<double> local_values(local_n);
    for (hsize_t i = 0; i < local_n; ++i)
        local_values[i] = static_cast<double>(offset + i);

    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
    assert(dxpl >= 0);
    assert(H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE) >= 0);

    assert(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
                    dxpl, local_values.data()) >= 0);
    assert(H5Pclose(dxpl) >= 0);
    assert(H5Sclose(memspace) >= 0);
    assert(H5Sclose(filespace) >= 0);
    assert(H5Dclose(dataset) >= 0);
    assert(H5Fclose(file) >= 0);

    MPI_Barrier(MPI_COMM_WORLD);

    fapl = H5Pcreate(H5P_FILE_ACCESS);
    assert(fapl >= 0);
    assert(H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL) >= 0);
    file = H5Fopen(filename, H5F_ACC_RDONLY, fapl);
    assert(file >= 0);
    H5Pclose(fapl);

    dataset = H5Dopen2(file, "values", H5P_DEFAULT);
    assert(dataset >= 0);
    filespace = H5Dget_space(dataset);
    assert(filespace >= 0);
    assert(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, nullptr,
                               &count, nullptr) >= 0);
    memspace = H5Screate_simple(1, &count, nullptr);
    assert(memspace >= 0);
    std::vector<double> read_values(local_n, -1.0);

    dxpl = H5Pcreate(H5P_DATASET_XFER);
    assert(dxpl >= 0);
    assert(H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE) >= 0);
    assert(H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
                   dxpl, read_values.data()) >= 0);

    for (hsize_t i = 0; i < local_n; ++i)
        assert(read_values[i] == static_cast<double>(offset + i));

    H5Pclose(dxpl);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dataset);
    H5Fclose(file);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        std::remove(filename);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}
