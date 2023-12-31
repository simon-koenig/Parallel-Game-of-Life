#include <iostream>
#include <mpi.h>
#include <cstddef>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    //*** Create cartesian topology for grids ***//
    int myrank, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    int dims[2] = {0, 0};
    const int ndims = 2;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Dims_create(numprocs, ndims, dims);
    MPI_Comm GRID_COMM;
    int bcs[2] = {1, 1};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, bcs, reorder, &GRID_COMM);
    int coords[2] = {};
    MPI_Barrier(GRID_COMM);
    MPI_Cart_coords(GRID_COMM, myrank, ndims, coords);

    int const n = 10;
    double matrix[n + 2][n + 2];
    const int t = 8;
    int sources[t], destinations[t];
    int target[] = {0, 1, 0, -1, -1, 0, 1, 0,
                    -1, 1, 1, 1, 1, -1, -1, -1};
    int vector[ndims];
    for (short i = 0; i < 8; i++)
    {
        MPI_Cart_coords(GRID_COMM, myrank, 2, vector);
        vector[0] += target[2 * i];
        vector[1] += target[2 * i + 1];
        MPI_Cart_rank(GRID_COMM, vector, &destinations[i]);
        MPI_Cart_coords(GRID_COMM, myrank, 2, vector);
        vector[0] -= target[2 * i];
        vector[1] -= target[2 * i + 1];
        MPI_Cart_rank(GRID_COMM, vector, &sources[i]);
    }
    MPI_Comm cartcomm;
    MPI_Dist_graph_create_adjacent(GRID_COMM,
                                   t, sources, MPI_UNWEIGHTED,
                                   t, destinations, MPI_UNWEIGHTED,
                                   MPI_INFO_NULL, reorder, &cartcomm);

    // Displacement Buffers
    MPI_Aint senddisp[8] = {1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Aint recvdisp[8] = {1, 1, 1, 1, 1, 1, 1, 1};
    // Type buffer
    MPI_Datatype sendtype[8];
    MPI_Datatype recvtype[8];
    // datatypes for matrix borders
    MPI_Datatype ROW, COL, COR;
    // Define the datatype for a column of the matrix
    MPI_Type_vector(n, 1, n, MPI_INT, &COL);
    MPI_Type_commit(&COL);
    MPI_Type_vector(n, 1, 1, MPI_INT, &ROW);
    MPI_Type_commit(&ROW);
    COR = MPI_INT;
    MPI_Type_commit(&COR);

    senddisp[0] = 1 * (n + 2) + n;
    sendtype[0] = COL;
    recvdisp[0] = 1 * (n + 2) + n + 1;
    recvtype[0] = COL;
    // ...
    senddisp[2] = 1 * (n + 2) + 1;
    sendtype[2] = ROW;
    recvdisp[2] = 1;
    recvtype[2] = ROW;
    // ...
    senddisp[4] = n * (n + 2) + 1;
    sendtype[4] = COR;
    recvdisp[4] = (n + 1) * (n + 2);
    recvtype[4] = COR;
    // ...
    // counts and byte offsets
    int sendcount[8];
    int recvcount[8];
    for (short i = 0; i < t; i++)
    {
        sendcount[i] = 1;
        recvcount[i] = 1;
        senddisp[i] *= sizeof(double);
        recvdisp[i] *= sizeof(double);
    }
    int iterate{1};
    for (short i = 0; i < 1; ++i)
    {
        // compute ...
        // update
        MPI_Neighbor_alltoallw(
            matrix, sendcount, senddisp, sendtype,
            matrix, recvcount, recvdisp, recvtype, cartcomm);
        int local = 1; // local convergence check
        MPI_Allreduce(&local, &iterate, 1, MPI_SHORT,
                      MPI_LAND, cartcomm);
    }
}
