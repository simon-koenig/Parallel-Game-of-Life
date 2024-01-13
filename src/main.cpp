#include <assert.h>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <array>
#include <cstring>
#ifdef USEMPI
#include <mpi.h>
#endif
#include "arguments.hpp"
#include "solver.hpp"

void StoreTimings(char* filename, int numproc, int repetitions, int resolution, int iterations, double* total, double* mean, double* max)
{   
    char datafile[120] = "./";
    strcat(datafile,filename);
    FILE *fp = fopen(datafile, "w");

    // store basic information about run
    fprintf(fp, "\n");
    fprintf(fp, "# processes = %i \n", numproc);
    fprintf(fp, "# repetitions = %i \n", repetitions);
    fprintf(fp, "# resolution = %i \n", resolution);
    fprintf(fp, "# iterations = %i \n", iterations);
    fprintf(fp, "\n");
    fprintf(fp, "summed time [s], mean time per processor [s], maximum time per processor [s] \n");    

    for(int i = 0; i < repetitions; i++)
    {
      // store values
      fprintf(fp, "%.6f,  %.6f, %.6f\n", total[i], mean[i], max[i]);
    }  
    
    fclose(fp);
}

int main(int argc, char *argv[])
{
  // parse command line arguments
  int rank = 0;
  int numproc = 1;
  auto repetitions = convertTo<int>(1, 10, argc, argv);
  auto resolution = convertTo<int>(2, 32, argc, argv);
  auto iterations = convertTo<int>(3, 1000, argc, argv);
  bool test_run = convertTo<bool>(4, 0, argc, argv);
  int ndims = 2;

  assert(resolution > 0);
  assert(iterations > 0);
  assert(numproc > 0);
  assert(repetitions > 0);

  double* total = (double*)malloc(sizeof(double)*repetitions);
  double* mean = (double*)malloc(sizeof(double)*repetitions);
  double* max = (double*)malloc(sizeof(double)*repetitions);

  std::array<double,3> timings;

  #ifdef USEMPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  if (rank == 0)
  {
    std::cout << "Starting Game of Life for " << iterations << " iterations with resolution of " << resolution << "living points." << std::endl;
  }


  for (int exp_counter = 1; exp_counter <= repetitions; exp_counter++)
  { 

    if (rank == 0)
    {
      std::cout << "|++++++++++++++++++++> " << "Starting Experiment " <<  exp_counter << "! |+++++++++++++++++++++|" << std::endl;
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

    timings = solve(resolution, iterations, rank, numproc, ndims, test_run);

    if (rank == 0)
    {
      total[exp_counter - 1] = timings[0];
      mean[exp_counter - 1] = timings[1];
      max[exp_counter - 1] = timings[2];
    }

    if (rank == 0)
    {
      std::cout << std::scientific << "|summed total runtime of all processes|= " << timings[0]<< " seconds" << std::endl;
      std::cout << std::scientific << "|maximum runtime of all processes|= " << timings[1] << " seconds" << std::endl;
      std::cout << std::scientific << "|average runtime per process|= " << timings[2] << " seconds per processor" << std::endl;
    };

    if (rank == 0)
    {    
      std::cout << "|++++++++++++++++++++> " << "Done with Experiment " << exp_counter << "!" << " |+++++++++++++++++++++|" << std::endl;
      fflush(stdout);
    };
  };

  if (rank == 0)
  {
    char filename[120];
    sprintf(filename, "%s%d%s%d%s%d%s%d%s", "./data/SENDRECV/P", numproc, "REPS", repetitions, "RES", resolution, "I", iterations, ".txt");
    StoreTimings(filename, numproc, repetitions, resolution, iterations, total, mean, max);
  }

  free(total);
  free(mean);
  free(max);

  #ifdef USEMPI
    MPI_Finalize();
  #endif

  return 0;
}
