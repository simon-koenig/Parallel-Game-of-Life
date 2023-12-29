#include <array>
#include <chrono>
#include <thread>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
// #ifdef USEMPI
#include <mpi.h>
// #endif
#include <assert.h>

// To change the number of processes in each direction change nx, ny
const int NX{10};
const int NY{10};

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

    //******** SENDING + RECEIVING VECTORS********//
    std::vector<int> NORTH_SEND(NX, 1);
    std::vector<int> NORTH_RECV(NX, 0);
    std::vector<int> SOUTH_SEND(NX, 2);
    std::vector<int> SOUTH_RECV(NX, 0);
    std::vector<int> EAST_SEND(NY, 3);
    std::vector<int> EAST_RECV(NY, 0);
    std::vector<int> WEST_SEND(NY, 4);
    std::vector<int> WEST_RECV(NY, 0);

    // Prepare data to be sent
    struct VectorStruct
    {
        int size;
        int *data;
    };

    const int data_size = 4;
    // Allocate space for the send array
    std::vector<VectorStruct> send_array(data_size);
    // Allocate space for the receive array
    std::vector<VectorStruct> recv_array(data_size);

    // Create an MPI datatype for VectorStruct
    MPI_Datatype vector_struct_type;
    MPI_Aint offsets[2] = {0, sizeof(int)};
    int block_lengths[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    MPI_Type_create_struct(2, block_lengths, offsets, types, &vector_struct_type);
    MPI_Type_commit(&vector_struct_type);

    // Initialize the data for each vector
    send_array[0].size = NORTH_SEND.size();
    send_array[0].data = NORTH_SEND.data();
    send_array[1].size = SOUTH_SEND.size();
    send_array[1].data = SOUTH_SEND.data();
    send_array[2].size = EAST_SEND.size();
    send_array[2].data = EAST_SEND.data();
    send_array[3].size = WEST_SEND.size();
    send_array[3].data = WEST_SEND.data();

    // Initialize the data for receive each vector
    recv_array[0].size = NORTH_RECV.size();
    recv_array[0].data = NORTH_RECV.data();
    recv_array[1].size = SOUTH_RECV.size();
    recv_array[1].data = SOUTH_RECV.data();
    recv_array[2].size = EAST_RECV.size();
    recv_array[2].data = EAST_RECV.data();
    recv_array[3].size = WEST_RECV.size();
    recv_array[3].data = WEST_RECV.data();
    // std::cout << "north send size " << send_array[0].size << std::endl;
    // std::cout << "north send data " << send_array[0].data[5] << std::endl;

    // Set up counts and displacements for MPI_Neighbor_alltoallw
    int sendcounts[data_size];
    int recvcounts[data_size];
    MPI_Aint sdispls[data_size];
    MPI_Aint rdispls[data_size];

    for (int i = 0; i < data_size; ++i)
    {
        sendcounts[i] = 1;
        recvcounts[i] = 1;
        sdispls[i] = i * sizeof(VectorStruct);
        rdispls[i] = i * sizeof(VectorStruct);
    };

    // Set up datatypes
    MPI_Datatype dtypes[4] = {vector_struct_type, vector_struct_type, vector_struct_type, vector_struct_type};

    // Rank buffers
    int north_rank;
    int south_rank;
    int east_rank;
    int west_rank;

    // // Get ranks of neighboring cells in the cartesian grid
    // MPI_Cart_shift(GRID_COMM, 0, 1, &east_rank, &west_rank);
    // MPI_Cart_shift(GRID_COMM, 1, 1, &south_rank, &north_rank);
    // int sources[4] = {north_rank, south_rank, east_rank, west_rank};
    // int destinations[4] = {north_rank, south_rank, east_rank, west_rank};
    // for (int i = 0; i < data_size; ++i)
    // {
    //     std::cout << "destinations" << destinations[i] << std::endl;
    //     std::cout << "sources" << sources[i] << std::endl;
    //     std::cout << "sdispls " << sdispls[i] << std::endl;
    // };
    // MPI_Comm GRAPH_COMM;
    // // Set up the graph topology
    // int indegree{4};  // north,south,east,west
    // int outdegree{4}; // north,south,east,west
    // // std::cout << "indegreee" << indegree << std::endl;
    // // std::cout << "outdegree" << outdegree << std::endl;
    // MPI_Barrier(GRID_COMM);

    // MPI_Dist_graph_create_adjacent(GRID_COMM, indegree, sources, MPI_UNWEIGHTED,
    //                                outdegree, destinations, MPI_UNWEIGHTED,
    //                                MPI_INFO_NULL, 0, &GRAPH_COMM);

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
    MPI_Comm GRAPH_COMM;
    MPI_Dist_graph_create_adjacent(GRID_COMM,
                                   t, sources, MPI_UNWEIGHTED,
                                   t, destinations, MPI_UNWEIGHTED,
                                   MPI_INFO_NULL, reorder, &GRAPH_COMM);

    // Print for debug
    for (int i = 0; i < data_size; ++i)
    {
        std::cout << "Rank " << myrank << " send: size =" << send_array[i].size
                  << ", data =";
        for (int j = 0; j < NX; ++j)
        {
            std::cout << send_array[i].data[j] << " ";
        }
        std::cout << std::endl;
    }
    fflush(stdout);
    MPI_Barrier(GRID_COMM); // Barrier before all send and recv
    // *************** Perform send *************** //
    MPI_Neighbor_alltoallw(&send_array, sendcounts, sdispls, dtypes,
                           &recv_array, recvcounts, rdispls, dtypes, GRAPH_COMM);
    MPI_Barrier(GRID_COMM);

    // std::cout << "Hi from after all to all " << std::endl;
    //  Print for debug
    // if (myrank == 0)
    // {
    //     for (int i = 0; i < data_size; ++i)
    //     {
    //         std::cout << "Print container received" << std::endl;
    //         std::cout << "Hi from struct: " << i << " of size:  " << recv_array[i].size << std::endl;
    //         for (int j = 0; j < NX; ++j)
    //         {
    //             std::cout << "Value " << recv_array[i].data[j] << "  " << std::endl;
    //         }
    //         std::cout << std::endl;
    //     };
    // }
    // fflush(stdout);

    // Print received data
    for (int i = 0; i < data_size; ++i)
    {
        std::cout << "Rank " << myrank << " received: size =" << recv_array[i].size
                  << ", data =";
        for (int j = 0; j < NX; ++j)
        {
            std::cout << recv_array[i].data[j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Hi from after printing " << std::endl;

    // fflush(stdout);

    MPI_Type_free(&vector_struct_type);
    MPI_Finalize();
    return 0;
}