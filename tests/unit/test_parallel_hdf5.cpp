#include <H5public.h>
#include <H5Fpublic.h>
#include <H5Ppublic.h>
#include <mpi.h>

#include <cassert>
#include <cstdio>
#include <string>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    assert(fapl >= 0);
    assert(H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL) >= 0);
    hid_t file = H5Fcreate("parallel_hdf5_capability_test.h5", H5F_ACC_TRUNC,
                           H5P_DEFAULT, fapl);
    assert(file >= 0);
    assert(H5Fclose(file) >= 0);
    assert(H5Pclose(fapl) >= 0);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        // The test only establishes MPI VFD availability and collective file
        // creation. Data distribution is deliberately left to the field I/O
        // integration work; this prevents a replicated write from masquerading
        // as distributed HDF5.
        std::remove("parallel_hdf5_capability_test.h5");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
