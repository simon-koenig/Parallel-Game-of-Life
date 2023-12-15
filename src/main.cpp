#include <assert.h>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#ifdef USEMPI
#include <mpi.h>
#endif
#include "arguments.hpp"
#include "solver.hpp"

int main(int argc, char *argv[])
{

  // parse command line arguments
  int rank = 0;
  auto numproc = convertTo<int>(1, 4, argc, argv);
  auto repetitions = convertTo<int>(2, 10, argc, argv);
  auto resolution = convertTo<int>(3, 32, argc, argv);
  auto iterations = convertTo<int>(4, 1000, argc, argv);
  int ndims = 2;

  assert(resolution > 0);
  assert(iterations > 0);
  assert(numproc > 0);
  assert(repetitions > 0);

for (int exp_counter = 1; exp_counter <= repetitions; exp_counter++)
{ 

  #ifdef USEMPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  if (rank == 0)
  {
    std::cout << "|++++++++++++++++++++> " << "Starting " << exp_counter << " Experiment!" << " |+++++++++++++++++++++|" << std::endl;
  };

  if (rank == 0)
  {
    std::cout << "numproc=" << numproc << std::endl;
    std::cout << "resolution=" << resolution << std::endl;
    std::cout << "iterations=" << iterations << std::endl;
  };

  if (rank == 0)
  {
    std::cout << "Start " << ndims << "-D Solver for [" << iterations << "] iterations with resolution of [" << resolution << "]" << std::endl;
  };

  std::array<double,3> timings;
  timings = solve(resolution, iterations, rank, numproc, ndims);

  if (rank == 0)
  {
    std::cout << std::scientific << "|summed total runtime of all processes|= " << timings[0]<< " seconds" << std::endl;
    std::cout << std::scientific << "|maximum runtime of all processes|= " << timings[1] << " seconds" << std::endl;
    std::cout << std::scientific << "|average runtime per process|= " << timings[2] << " seconds per processor" << std::endl;
  }

  if (rank == 0)
  {    
    std::cout << "|++++++++++++++++++++> " << "Done with Experiment " << exp_counter << "!" << " |+++++++++++++++++++++|" << std::endl;
    fflush(stdout);
  };
  #ifdef USEMPI
    MPI_Finalize();
  #endif
};

  return 0;
}
