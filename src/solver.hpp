#pragma once
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
#include "game_sequential.hpp"

template <typename Type>
class MatrixView
{
private:
  std::vector<Type> &v;
  MatrixView(const MatrixView &);
  MatrixView &operator=(const MatrixView &);

public:
  const int N, M;
  MatrixView(std::vector<Type> &v, int N, int M) : v(v), N(N), M(M)
  {
    assert(static_cast<int>(v.size()) / N == M);
  }
  Type &set(int i, int j) { return v[i + N * j]; }
  const Type &get(int i, int j) { return v[i + N * j]; }
  Type &set(int n) { return v[n]; }
  const Type &get(int n) { return v[n]; }

  // function to access the underlying vector
  std::vector<Type> &getVector()
  {
    return v;
  }
  const std::vector<Type> &getVector() const
  {
    return v;
  }
};


template <typename Type>
void printMatrices(int myrank, int mpi_numproc, int NY, int NX, MatrixView<Type> &solutionView, bool excl_ghost_layer=true)
{
  int offset = 1;
  if (excl_ghost_layer == false){
    offset = 0;
  }

  // sleep to make sure that everything gets printed correctly
  this_thread::sleep_for(chrono::milliseconds(100 * myrank));

  for (int i = 0; i < mpi_numproc; i++)
  {
    if (myrank == i)
    {
      printf("%d: matrix\n", myrank);
      for (int j = offset; j != NY - offset; ++j) // Indexing to onyl print system matrix w/o ghost layers
      {
        for (int i = offset; i != NX - offset; ++i)
        {
          //std::cout << solutionView.get(i, j) << " ";
          std::cout << (solutionView.get(i, j) == 1 ? "■" : "□") << " ";
        }
        std::cout << ("\n");
      }
    }
    fflush(stdout);
  }
}

template <typename Type>
std::vector<int> getGrid(int myrank, int numprocs, int (&dims)[2], int (&coords)[2], int NX, int NY, MatrixView<Type> &solutionView, MPI_Comm COMM)
{
  // collect submatrices to assemble full matrix for correctness testing
  int full_matrix_size = dims[0]*dims[1]*(NX-2)*(NY-2);

  std::vector<int> submatrix_collection(NX * NY * dims[0] * dims[1]);
  //std::cout << "SIZES: " << submatrix_collection.size() << " " << solutionView.getVector().size() << std::endl;

  // cheap way to "make sure" that gather is called in the correct order
  this_thread::sleep_for(chrono::milliseconds(100 * myrank));
  
  // get coordinates of submatrices
  std::vector<int> all_coords(2 * numprocs);
  MPI_Gather(&coords, dims[0] , MPI_INT, 
             all_coords.data() , dims[0] , MPI_INT , 0 , COMM);

  // get submatrices
  MPI_Gather(solutionView.getVector().data(), solutionView.getVector().size(), MPI_INT,
             submatrix_collection.data(), solutionView.getVector().size(), MPI_INT,
             0, COMM);
  
  MPI_Barrier(COMM);

  std::vector<int> fm(full_matrix_size);
  MatrixView<int> full_matrix (fm, dims[0]*(NX-2), dims[1]*(NY-2));

  if (myrank == 0) {
    /*
    std::cout << "All Coordinates:\n";
    for (int i = 0; i < numprocs; ++i) {
      std::cout << "Process " << i << ": (" << all_coords[2 * i] << ", " << all_coords[2 * i + 1] << ")\n";
    }*/

    int count = 1;
    for (int p = 0; p < dims[0]*dims[1]; p++){
      //std::cout << all_coords[2 * p] << all_coords[2 * p + 1] << std::endl;
      count += NX;
      int y = 0;
      for (int j = 1 + p*(NY-2); j < 1 + (p + 1)*(NY-2); j++){
        int x = 0;
        for (int i = j*(NX) + 1; i < (j + 1)*(NX) - 1; i++){
          full_matrix.set(x+(NX-2)*all_coords[2 * p], y+(NY-2)*all_coords[2 * p + 1]) = submatrix_collection[count];
          x++;
          count++;
        }
        y++;
        count += 2;
      }
      count += NX;
    }
    //printMatrices(myrank, numprocs, dims[0]*(NX-2), dims[1]*(NY-2), full_matrix, false);
  } 
  return full_matrix.getVector();
}

// Cell Naming for easy access
enum Cell
{
  UNKNOWN = 0,
  NORTH = 1,
  SOUTH = 2,
  EAST = 3,
  WEST = 4,
  NORTHEAST = 5,
  SOUTHEAST = 6,
  SOUTHWEST = 7,
  NORTHWEST = 8

};

// Function to get seeded random value
uint8_t get_random_value(const int row_id, const int col_id, const int n,
                         const int seed)
{
  uint8_t r = 0;
  int my_seed = seed + row_id * n + col_id;
  // printf("my_seed: %d\n", my_seed);
  srand(my_seed);
  // rand();
  // printf("rand=%d\n", rand());
  rand();
  r = rand() % 2;
  return r;
}

std::array<double, 3> solve(int resolution, int iterations, int mpi_rank,
                            int mpi_numproc, int ndims)
{

  // Solver Template stops here, MPI subdomain solve starts
  //  #ifdef USEMPI

  int myrank, numprocs;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  // Build Cartesian Grid
  //  std::cout << "Init Processor Rank " << myrank << std::endl;
  int dims[2] = {0, 0};

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Dims_create(numprocs, ndims, dims);
  MPI_Comm GRID_COMM;
  int bcs[2] = {1, 1};
  int reorder = 1;
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, bcs, reorder, &GRID_COMM);
  int coords[2] = {};
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Cart_coords(GRID_COMM, myrank, ndims, coords);

  int width_x, width_y, last_width_x, last_width_y;
  int modulo_x, modulo_y;
  int NY = 0;
  int NX = 0;
  bool is_subdomain = coords[0] + 1 != dims[0] && coords[1] + 1 != dims[1];
  bool is_x_last_subdomain = coords[0] + 1 == dims[0] && coords[1] + 1 != dims[1];
  bool is_y_last_subdomain = coords[0] + 1 != dims[0] && coords[1] + 1 == dims[1];
  bool is_xy_last_subdomain = coords[0] + 1 == dims[0] && coords[1] + 1 == dims[1];

  width_x = (resolution) / dims[0];
  width_y = (resolution) / dims[1];

  modulo_x = (resolution) % dims[0];
  modulo_y = (resolution) % dims[1];

  last_width_x = width_x + modulo_x;
  last_width_y = width_y + modulo_y;

  if (myrank == 0)
  {
    std::cout << "numprocs = " << numprocs << std::endl;
    std::cout << "dims=(" << dims[0] << "," << dims[1] << ")" << std::endl;
  }

  if (is_subdomain)
  {
    NX = width_x + 2;
    NY = width_y + 2;
    std::cout << "NX: " << NX << "  NY:   " << NY << std::endl;
  }
  else if (is_x_last_subdomain)
  {
    NX = last_width_x + 2;
    NY = width_y + 2;
  }
  else if (is_y_last_subdomain)
  {
    NX = width_x + 2;
    NY = last_width_y + 2;
  }
  else if (is_xy_last_subdomain)
  {
    NX = last_width_x + 2;
    NY = last_width_y + 2;
  }

  std::vector<int> domain(NX * NY, Cell::UNKNOWN);
  MatrixView<int> domainView(domain, NX, NY);

  // TODO: Adapt this, there are no boundary conditions. Our domain is periodical in alle directions
  for (int i = 1; i != NX - 1; ++i)
  {
    domainView.set(i, 0) = Cell::NORTH;

    domainView.set(i, NY - 1) = Cell::SOUTH;
  }

  for (int j = 1; j != NY - 1; ++j)
  {
    domainView.set(NX - 1, j) = Cell::EAST;
    domainView.set(0, j) = Cell::WEST;
  }

  // Set Diagonal Elements
  domainView.set(0, 0) = Cell::NORTHWEST;
  domainView.set(NX - 1, 0) = Cell::NORTHEAST;
  domainView.set(0, NY - 1) = Cell::SOUTHWEST;
  domainView.set(NX - 1, NY - 1) = Cell::SOUTHEAST;

  //  std::cout << "Init Done Processor Rank " << myrank << std::endl;

  MPI_Barrier(GRID_COMM);
  //	std::cout << "coords: " << coords[0] << "," << coords[1] << " rank: " << myrank << " size: " << NX << "," << NY << std::endl;

  // TODO: Play the game of life here. We only have cross stencil (4 neighbours), but need star stencil (8 neighbours)
  auto game = [](std::vector<int> &sol, std::vector<int> &sol2,
                 int NX, int NY)
  {
    MatrixView<int> solView(sol, NX, NY);
    MatrixView<int> sol2View(sol2, NX, NY);
    // MatrixView<int> rhsView(rhs, NX, NY);

    for (int i = 1; i != NX - 1; ++i)
    {
      for (int j = 1; j != NY - 1; ++j)
      {
        // Get neighbor values
        int center = solView.get(i, j);
        int north = solView.get(i, j - 1);
        int south = solView.get(i, j + 1);
        int east = solView.get(i + 1, j);
        int west = solView.get(i - 1, j);
        int northeast = solView.get(i + 1, j - 1);
        int southeast = solView.get(i + 1, j + 1);
        int southwest = solView.get(i - 1, j + 1);
        int northwest = solView.get(i - 1, j - 1);
        int sum{0};
        sum = north + south + east + west +
              southeast + southwest + northeast + northwest;
        // std::cout << "Hey from (i,j) = ( " << i << " , " << j << " ) , center = " << center
        //           << " north = " << north << " south = " << south << " east = " << east << " west  =  " << west << std::endl;
        // std::cout << "Hey from (i,j) = ( " << i << " , " << j << " ) , sum = " << sum << std::endl;

        // (huh huh huh huh) Staying alive condition
        if (center == 1 && sum == 2)
        {
          sol2View.set(i, j) = 1;
        }
        // (huh huh huh huh) Staying alive condition
        else if (center == 1 && sum == 3)
        {
          sol2View.set(i, j) = 1;
        }
        // Reborn condition
        else if (center == 0 && sum == 3)
        {
          sol2View.set(i, j) = 1;
        }
        // Death condition
        else
        {
          sol2View.set(i, j) = 0;
        }
      }
    }
    sol.swap(sol2);
  };

  // TODO: Build initial grid, set randomly to 0 or 1
  //  solution approximation starting with boundary initialized to dirichlet
  //  conditions, else 0

  // Buffer for own coordinates to get the rank of diagonal coords
  int own_coords[2];
  int diag_coords[2];

  //
  // Initialise Grid
  //

  int m_offset_r = own_coords[1] * NY;
  int m_offset_c = own_coords[0] * NX;
  std::vector<int> solution(NX * NY, 0);
  MatrixView<int> solutionView(solution, NX, NY);
  for (int j = 1; j != NY - 1; ++j)
  {
    for (int i = 1; i != NX - 1; ++i)
    {
      solutionView.set(i, j) = get_random_value(m_offset_r + i, m_offset_c + j, 2, 10);
    }
  };

  // Init Send and Recv Buffers
  std::vector<int> solution2 = solution;
  std::vector<int> NORTH_SEND(NX, 0);
  std::vector<int> NORTH_RECV(NX, 0);
  std::vector<int> SOUTH_SEND(NX, 0);
  std::vector<int> SOUTH_RECV(NX, 0);
  std::vector<int> EAST_SEND(NY, 0);
  std::vector<int> EAST_RECV(NY, 0);
  std::vector<int> WEST_SEND(NY, 0);
  std::vector<int> WEST_RECV(NY, 0);
  int NORTHEAST_RECV{0};
  int NORTHEAST_SEND{0};
  int NORTHWEST_RECV{0};
  int NORTHWEST_SEND{0};
  int SOUTHWEST_RECV{0};
  int SOUTHWEST_SEND{0};
  int SOUTHEAST_RECV{0};
  int SOUTHEAST_SEND{0};

  // Rank buffers
  int north_rank;
  int south_rank;
  int east_rank;
  int west_rank;
  int north_east_rank;
  int south_east_rank;
  int north_west_rank;
  int south_west_rank;

  // Get ranks of neighboring cells in the cartesian grid
  MPI_Cart_shift(GRID_COMM, 0, 1, &east_rank, &west_rank);
  MPI_Cart_shift(GRID_COMM, 1, 1, &south_rank, &north_rank);

  MPI_Cart_coords(GRID_COMM, myrank, 2, own_coords); // Get own coords

  // North East Diagonal
  diag_coords[0] = (own_coords[0] - 1 + dims[0]) % dims[0];
  diag_coords[1] = (own_coords[1] + 1) % dims[1];
  MPI_Cart_rank(GRID_COMM, diag_coords, &north_east_rank); // Norht east rank now holds the rank processor of the processor north east

  // South East Diagonal
  diag_coords[0] = (own_coords[0] + 1) % dims[0];
  diag_coords[1] = (own_coords[1] + 1) % dims[1];
  MPI_Cart_rank(GRID_COMM, diag_coords, &south_east_rank); // South east rank now holds the rank processor of the processor north east

  // North west Diagonal
  diag_coords[0] = (own_coords[0] - 1 + dims[0]) % dims[0];
  diag_coords[1] = (own_coords[1] - 1 + dims[1]) % dims[1];
  MPI_Cart_rank(GRID_COMM, diag_coords, &north_west_rank); // North west rank now holds the rank processor of the processor north east

  // South west Diagonal
  diag_coords[0] = (own_coords[0] + 1) % dims[0];
  diag_coords[1] = (own_coords[1] - 1 + dims[1]) % dims[1];
  MPI_Cart_rank(GRID_COMM, diag_coords, &south_west_rank); // South west rank now holds the rank processor of the processor north east

  // get whole grid before parallel run
  //int grid_size = dims[0]*dims[1]*(NX-2)*(NY-2);
  //std::vector<int> grid_before(grid_size);
  //grid_before = getGrid(myrank, numprocs, dims, coords, NX, NY, solutionView, GRID_COMM);

  MPI_Barrier(GRID_COMM);

  //if (myrank == 0)
  //{
  //std::cout << "\nSolve Game of Life using 8 point stencil:" << std::endl
  //            << std::endl;
  //std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
  //  std::cout << "++++++ Before Lifetime +++++++++++++++" << std::endl;
  //  std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
  //  fflush(stdout);
  //}

  // Print Grid before lifetime
  // printMatrices(myrank, mpi_numproc, NY, NX, solutionView);
  // Wait for all processes to be finished
  MPI_Barrier(GRID_COMM);

  // START THE GAME OF LIFE
  auto start = MPI_Wtime();
  
  for (int iter = 0; iter <= iterations; ++iter)
  {

    // Preparing ghost layer for sending
    // => Handling ghost layers the data to be sent

    for (int j = 0; j != NY; ++j)
    {
      for (int i = 0; i != NX; ++i)
      {
        if (domainView.get(i, j) == Cell::NORTH)
        {                                             //&& dims[1] != 1) {
          NORTH_SEND[i] = solutionView.get(i, j + 1); // x
        }
        else if (domainView.get(i, j) == Cell::SOUTH)
        {                                             //&& dims[1] != 1){
          SOUTH_SEND[i] = solutionView.get(i, j - 1); // x
        }
        else if (domainView.get(i, j) == Cell::EAST)
        {
          EAST_SEND[j] = solutionView.get(i - 1, j);
        }
        else if (domainView.get(i, j) == Cell::WEST)
        {
          WEST_SEND[j] = solutionView.get(i + 1, j);
        };
      }
    }

    if (ndims != 1) // If grid is only one dimensional
    {
      MPI_Sendrecv(NORTH_SEND.data(), NX, MPI_INT, north_rank, 1,
                   SOUTH_RECV.data(), NX, MPI_INT, south_rank, 1, GRID_COMM, MPI_STATUS_IGNORE);
      MPI_Barrier(GRID_COMM);
      MPI_Sendrecv(SOUTH_SEND.data(), NX, MPI_INT, south_rank, 2,
                   NORTH_RECV.data(), NX, MPI_INT, north_rank, 2, GRID_COMM, MPI_STATUS_IGNORE);
      MPI_Barrier(GRID_COMM);
    };
    MPI_Sendrecv(WEST_SEND.data(), NY, MPI_INT, west_rank, 3,
                 EAST_RECV.data(), NY, MPI_INT, east_rank, 3, GRID_COMM, MPI_STATUS_IGNORE);
    MPI_Barrier(GRID_COMM);
    MPI_Sendrecv(EAST_SEND.data(), NY, MPI_INT, east_rank, 4,
                 WEST_RECV.data(), NY, MPI_INT, west_rank, 4, GRID_COMM, MPI_STATUS_IGNORE);
    MPI_Barrier(GRID_COMM);

    // TODO: DIAGONAL SEND
    // int NORTHEAST_RECV{0};

    // int NORTHWEST_RECV{0};
    // int NORTHWEST_SEND{0};

    // int SOUTHWEST_SEND{0};
    // int SOUTHEAST_RECV{0};
    // int SOUTHEAST_SEND{0};

    // Northeast Send, Southwest Recv
    NORTHEAST_SEND = solutionView.get(NX - 1, NY - 1);

    MPI_Sendrecv(&NORTHEAST_SEND, 1, MPI_INT, north_east_rank, 5,
                 &SOUTHWEST_RECV, 1, MPI_INT, south_west_rank, 5, GRID_COMM, MPI_STATUS_IGNORE);
    solutionView.set(0, 0) = SOUTHWEST_RECV;
    MPI_Barrier(GRID_COMM);

    // TODO: Northwest Send, Southeast Recv
    NORTHWEST_SEND = solutionView.get(0, NY - 1);

    MPI_Sendrecv(&NORTHWEST_SEND, 1, MPI_INT, north_east_rank, 6,
                 &SOUTHEAST_RECV, 1, MPI_INT, south_west_rank, 6, GRID_COMM, MPI_STATUS_IGNORE);
    solutionView.set(NX - 1, 0) = SOUTHEAST_RECV;
    MPI_Barrier(GRID_COMM);

    // TODO: Southeast Send, Northwest Recv
    SOUTHEAST_SEND = solutionView.get(NX - 1, 0);

    MPI_Sendrecv(&SOUTHEAST_SEND, 1, MPI_INT, north_east_rank, 6,
                 &NORTHWEST_RECV, 1, MPI_INT, south_west_rank, 6, GRID_COMM, MPI_STATUS_IGNORE);
    solutionView.set(0, NY - 1) = NORTHWEST_RECV;
    MPI_Barrier(GRID_COMM);

    // TODO: Southwest Send, Northeast Recv

    SOUTHWEST_SEND = solutionView.get(0, 0);

    MPI_Sendrecv(&SOUTHWEST_SEND, 1, MPI_INT, north_east_rank, 6,
                 &NORTHEAST_RECV, 1, MPI_INT, south_west_rank, 6, GRID_COMM, MPI_STATUS_IGNORE);
    solutionView.set(NX - 1, NY - 1) = NORTHEAST_RECV;
    MPI_Barrier(GRID_COMM);

    // Update domain data with data received by the ghost layers
    //
    // Non Diagonal Send
    for (int j = 0; j != NY; ++j)
    {
      for (int i = 0; i != NX; ++i)
      {
        if (domainView.get(i, j) == Cell::NORTH)
        { //&& dims[1] != 1){
          solutionView.set(i, j) = NORTH_RECV[i];
        }
        else if (domainView.get(i, j) == Cell::SOUTH)
        { //&& dims[1] != 1){
          solutionView.set(i, j) = SOUTH_RECV[i];
        }
        else if (domainView.get(i, j) == Cell::EAST)
        {
          solutionView.set(i, j) = EAST_RECV[j];
        }
        else if (domainView.get(i, j) == Cell::WEST)
        {
          solutionView.set(i, j) = WEST_RECV[j];
        };
      }
    }

    // Diagonal Send

    // Northeast Send, Southwest Recv
    NORTHEAST_SEND = solutionView.get(NX - 2, 1);

    MPI_Sendrecv(&NORTHEAST_SEND, 1, MPI_INT, north_east_rank, 5,
                 &SOUTHWEST_RECV, 1, MPI_INT, south_west_rank, 5, GRID_COMM, MPI_STATUS_IGNORE);
    // Southwest Recv
    solutionView.set(0, NY - 1) = SOUTHWEST_RECV;
    MPI_Barrier(GRID_COMM);

    // TODO: Northwest Send
    NORTHWEST_SEND = solutionView.get(1, 1);

    MPI_Sendrecv(&NORTHWEST_SEND, 1, MPI_INT, north_east_rank, 6,
                 &SOUTHEAST_RECV, 1, MPI_INT, south_west_rank, 6, GRID_COMM, MPI_STATUS_IGNORE);
    // Southeast Recv
    solutionView.set(NX - 1, NY - 1) = SOUTHEAST_RECV;
    MPI_Barrier(GRID_COMM);

    // TODO: Southeast Send
    SOUTHEAST_SEND = solutionView.get(NX - 2, NY - 2);

    MPI_Sendrecv(&SOUTHEAST_SEND, 1, MPI_INT, north_east_rank, 6,
                 &NORTHWEST_RECV, 1, MPI_INT, south_west_rank, 6, GRID_COMM, MPI_STATUS_IGNORE);
    // Northwest Recv
    solutionView.set(0, 0) = NORTHWEST_RECV;
    MPI_Barrier(GRID_COMM);

    // TODO: Southwest Send
    SOUTHWEST_SEND = solutionView.get(1, NY - 2);

    MPI_Sendrecv(&SOUTHWEST_SEND, 1, MPI_INT, north_east_rank, 6,
                 &NORTHEAST_RECV, 1, MPI_INT, south_west_rank, 6, GRID_COMM, MPI_STATUS_IGNORE);
    // Northeast Recv
    solutionView.set(NX - 1, 0) = NORTHEAST_RECV;

    MPI_Barrier(GRID_COMM);

    // Update own domain
    game(solution, solution2, NX, NY);
  };
  // Wait for all processes to be finished
  MPI_Barrier(GRID_COMM);

  {
    auto stop = MPI_Wtime();
    auto seconds = stop - start;

    double seconds_sum;
    double seconds_max;

    MPI_Allreduce(&seconds, &seconds_sum, 1, MPI_DOUBLE, MPI_SUM, GRID_COMM);
    MPI_Allreduce(&seconds, &seconds_max, 1, MPI_DOUBLE, MPI_MAX, GRID_COMM);

    MPI_Barrier(GRID_COMM);
    //
    // Print results
    //

    //if (myrank == 0)
    //{
    //  std::cout << "+++++++++++++++++++++++++++++++++++" << std::endl;
    //  std::cout << "++++++After Lifetime+++++++++++++++" << std::endl;
    //  std::cout << "+++++++++++++++++++++++++++++++++++" << std::endl;
    //  fflush(stdout);
    //}

    // printMatrices(myrank, mpi_numproc, NY, NX, solutionView);

    // whole grid 
    //std::vector<int> grid_after(grid_size);
    //grid_after = getGrid(myrank, numprocs, dims, coords, NX, NY, solutionView, GRID_COMM);
      
    // sequential run
    //if (myrank == 0){
    //  std::cout << "+++++++++++++++++ INITIAL MATRIX ++++++++++++++++++++\n" << std::endl;
    //  printGrid(grid_before, resolution, resolution);
    //  std::cout << "++++++++++++++++ AFTER PARALLEL RUN +++++++++++++++++\n" << std::endl;
    //  printGrid(grid_after, resolution, resolution);
    //  std::cout << "+++++++++++++++ AFTER SEQUENTIAL RUN ++++++++++++++++\n" << std::endl;
    // run_sequential(grid_before, resolution, resolution, iterations);
    //}

    std::array<double, 3> timings;
    timings[0] = seconds_sum;
    timings[1] = seconds_max;
    timings[2] = seconds_sum / numprocs;
    return timings;
  };
}
