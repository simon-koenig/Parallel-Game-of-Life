#include <mpi.h>
#include <iostream>

// Simple mpi test program
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << "Hello from process " << rank << " out of " << size << " processes." << std::endl;

    MPI_Finalize();
    return 0;
}