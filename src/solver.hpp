#pragma once
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
// #ifdef USEMPI
#include <mpi.h>
// #endif
#include <assert.h>

template <typename Type>
class MatrixView
{
private:
  std::vector<Type> &v;
  MatrixView(const MatrixView &);
  MatrixView &operator=(const MatrixView &);

public:
  const size_t N, M;
  MatrixView(std::vector<Type> &v, size_t N, size_t M) : v(v), N(N), M(M)
  {
    assert(v.size() / N == M);
  }
  Type &set(size_t i, size_t j) { return v[i + N * j]; }
  const Type &get(size_t i, size_t j) { return v[i + N * j]; }
  Type &set(size_t n) { return v[n]; }
  const Type &get(size_t n) { return v[n]; }
};

template <typename Type>
void printMatrices(int myrank, int mpi_numproc, int NY, int NX, MatrixView<Type> &solutionView)
{
  for (int i = 0; i < mpi_numproc; i++)
  {
    if (myrank == i)
    {
      printf("%d: matrix\n", myrank);
      for (size_t j = 0; j != NY; ++j)
      {
        for (size_t i = 0; i != NX; ++i)
        {
          std::cout << solutionView.get(i, j);
        }
        std::cout << ("\n");
      }
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

enum Cell
{
  UNKNOWN = 0,
  DIR = 1,
  SOUTH = 2,
  NORTH = 3,
  EAST = 4,
  WEST = 5
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

void solve(size_t resolution, size_t iterations, int mpi_rank,
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
  if (ndims == 1)
  {
    dims[1] = 1;
  }
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
  size_t NY = 0;
  size_t NX = 0;
  bool is_subdomain = coords[0] + 1 != dims[0] && coords[1] + 1 != dims[1];
  bool is_x_last_subdomain = coords[0] + 1 == dims[0] && coords[1] + 1 != dims[1];
  bool is_y_last_subdomain = coords[0] + 1 != dims[0] && coords[1] + 1 == dims[1];
  bool is_xy_last_subdomain = coords[0] + 1 == dims[0] && coords[1] + 1 == dims[1];

  width_x = (resolution - 2) / dims[0];
  width_y = (resolution - 2) / dims[1];

  modulo_x = (resolution - 2) % dims[0];
  modulo_y = (resolution - 2) % dims[1];

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
  for (size_t i = 1; i != NX - 1; ++i)
  {

    domainView.set(i, 0) = Cell::SOUTH;

    domainView.set(i, NY - 1) = Cell::NORTH;
  }

  for (size_t j = 1; j != NY - 1; ++j)
  {
    domainView.set(0, j) = Cell::EAST;
    domainView.set(NX - 1, j) = Cell::WEST;
  }

  //  std::cout << "Init Done Processor Rank " << myrank << std::endl;

  MPI_Barrier(GRID_COMM);
  //	std::cout << "coords: " << coords[0] << "," << coords[1] << " rank: " << myrank << " size: " << NX << "," << NY << std::endl;

  // TODO: Play the game of life here. We only have cross stencil (4 neighbours), but need star stencil (8 neighbours)
  auto SolverJacobi = [](std::vector<int> &sol, std::vector<int> &sol2,
                         size_t NX, size_t NY)
  {
    MatrixView<int> solView(sol, NX, NY);
    MatrixView<int> sol2View(sol2, NX, NY);
    // MatrixView<int> rhsView(rhs, NX, NY);

    for (size_t j = 1; j != NY - 1; ++j)
    {
      for (size_t i = 1; i != NX - 1; ++i)
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
        int sum = north + south + east + west +
                  southeast + southwest + northeast + northwest;

        // (huh huh huh huh) Staying alive condition
        if ((center == 1 && sum == 2) || (center == 1 && sum == 3))
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

  std::vector<int> solution(NX * NY, 0);
  MatrixView<int> solutionView(solution, NX, NY);
  for (size_t j = 0; j != NY; ++j)
  {
    for (size_t i = 0; i != NX; ++i)
    {
      solutionView.set(i, j) = get_random_value(myrank + j, myrank + i, ndims, 10);
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

  // Buffer for own coordinates to get the rank of diagonal coords
  int own_coords[2];
  int diag_coords[2];
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

  if (myrank == 0)
  {
    std::cout << "solve Game of Life using 8 point stencil" << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "++++++ Before Lifetime +++++++++++++++" << std::endl;
    std::cout << "++++++++++++++++++++++++++++++++++++++" << std::endl;
    fflush(stdout);
  }

  // Print Grid before lifetime
  printMatrices(myrank, mpi_numproc, NY, NX, solutionView);

  auto start = std::chrono::high_resolution_clock::now();
  for (size_t iter = 0; iter <= iterations; ++iter)
  {

    // Update own domain
    SolverJacobi(solution, solution2, NX, NY);

    // Preparing ghost layer for sending
    // => Handling ghost layers the data to be sent

    for (size_t j = 0; j != NY; ++j)
    {
      for (size_t i = 0; i != NX; ++i)
      {
        if (domainView.get(i, j) == Cell::NORTH)
        {                                             //&& dims[1] != 1) {
          NORTH_SEND[i] = solutionView.get(i, j - 1); // x
        }
        else if (domainView.get(i, j) == Cell::SOUTH)
        {                                             //&& dims[1] != 1){
          SOUTH_SEND[i] = solutionView.get(i, j + 1); // x
        }
        else if (domainView.get(i, j) == Cell::EAST)
        {
          EAST_SEND[j] = solutionView.get(i + 1, j);
        }
        else if (domainView.get(i, j) == Cell::WEST)
        {
          WEST_SEND[j] = solutionView.get(i - 1, j);
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

    // TODO: DIAGONAL SEND
    // Northeast Send, Southwest Recv
    MPI_Sendrecv(solutionView.get(NX - 1, NY - 1), 1, MPI_INT, north_east_rank, 5,
                 solutionView.set(0, 0), 1, MPI_INT, south_west_rank, 5, GRID_COMM, MPI_STATUS_IGNORE);
    MPI_Barrier(GRID_COMM);

    MPI_Barrier(GRID_COMM); // All processes have to wait for the others after each iteration for solution to be valid

    // Update domain data with data received by the ghost layers
    for (size_t j = 0; j != NY; ++j)
    {
      for (size_t i = 0; i != NX; ++i)
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
  };
  // Wait for all processes to be finished
  MPI_Barrier(GRID_COMM);

  {
    auto stop = std::chrono::high_resolution_clock::now();
    auto seconds = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();

    double seconds_sum;

    MPI_Reduce(&seconds, &seconds_sum, 1, MPI_DOUBLE, MPI_SUM, 0, GRID_COMM);

    if (myrank == 0)
    {
      std::cout << std::scientific << "|total runtime|= " << seconds_sum << " seconds" << std::endl;
      std::cout << std::scientific << "|average runtime per process|= " << seconds / numprocs << " seconds per processor" << std::endl;
    }
    //  #endif

    //
    // Print results
    //

    if (myrank == 0)
    {
      std::cout << "+++++++++++++++++++++++++++++++++++" << std::endl;
      std::cout << "++++++After Lifetime+++++++++++++++" << std::endl;
      std::cout << "+++++++++++++++++++++++++++++++++++" << std::endl;
      fflush(stdout);
    }

    printMatrices(myrank, mpi_numproc, NY, NX, solutionView);
  };
}
