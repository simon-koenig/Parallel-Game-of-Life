Project about Parallel-Game-of-Life for "High Performance Computing", winter term 2023 at TU WIEN.

Contributors:
* Simon König, 11702826
* Viktor Beck, 11713110
* Michael Ketter, 11701297

# Project descripton
The Game of Life is a simple time-stepping algorithm, which updates a binary 2D matrix according to predefined rules. 
Those rules usually concern only the value of the updated cell and its direct neighbours, making the algorithm local per definition.
Due to this locality, the algorithm is trivially parallelizable by splitting the 2D domain into submatrices.
Each submatrice is assigned to an available processor, which solely performs the updates of its submatrix.
Between the update steps, each processor has to communicate its border values to neighbouring processors. 
However, this communication step creates a bottleneck for parallel speed-up when increasing the processor number. 

The aim of this project is to implement a sequential and multiple parallel versions (utilizing the features in OpenMPI) of the game of life.
To assess the versions performance, the code should be benchmarked in a weak and strong scaling parallel set-up.
Findings, explanations and possible improvements should be stated and discussed in a short report.

# Dependencies
To run all the plotting scripts one needs to have the following packages installed (in python3):
* argparse
* os
* pathlib import Path
* itertools
* matplotlib
* numpy

Additionally, mpic++ is needed for compiling the C++-code into working binaries.

# Usage
## Running jobs locally 
The supplied makefile does not support running the code locally (without having slurm set-up).
However, with following commands the code can be compiled and run even locally.
Running the parallel version utilizing only ```MPI_Sendrecv```: 
```sh
mpic++ ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI -std=c++14 -O3 -Wall -pedantic -march=native
mpirun -n NPROCS --use-hwthread-cpus ./bin/GameOfLifeMPI REPETITION RESOLUTION ITERATIONS CORRECTNESS_TEST
```

Running the parallel version utilizing ```MPI_Neighbor_alltoallw```: 
```sh
mpic++ ./src/mainATA.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI -std=c++14 -O3 -Wall -pedantic -march=native
mpirun -n NPROCS --use-hwthread-cpus ./bin/GameOfLifeMPI REPETITION RESOLUTION ITERATIONS CORRECTNESS_TEST
```

## Running jobs on a cluster (e.g. hydra):
The makefile supports running the code on a cluster where job managment is done via slurm.
However, to import the MPI library one has to still run
```sh
spack load openmpi@4.1.5
```
before starting a benchmark.
To recreate the benchmark results shown in the report, one has to call
```sh
make StrongScalingExperiment IMPLEMENTATION=X
```
respectively,
```sh
make WeakScalingExperiment IMPLEMENTATION=X
```

After generating the benchmark data, one can generate the result plots by running:
```sh
make EvaluateStrongScalingExperiment IMPLEMENTATION=X
```
respectively,
```sh
make EvaluateWeakScalingExperiment IMPLEMENTATION=X
```

For further information about variation options, please refer to comments in the makefile and code.

# Folder structure
 AssignmentDescription.pdf
├── LICENSE
├── Makefile
├── README.md
├── bin
│   └── GameOfLifeMPI
├── data
│   ├── ALLTOALL
│   │   ├── P1024REPS20RES10240I50.txt
│   │   ├── P1024REPS20RES1024I50.txt
│   │   ├── P1024REPS20RES32768I50.txt
│   │   ├── P1REPS20RES10240I50.txt
│   │   ├── P1REPS20RES1024I50.txt
│   │   ├── P256REPS20RES10240I50.txt
│   │   ├── P256REPS20RES1024I50.txt
│   │   ├── P256REPS20RES16384I50.txt
│   │   ├── P32REPS20RES10240I50.txt
│   │   ├── P32REPS20RES1024I50.txt
│   │   ├── P32REPS20RES5793I50.txt
│   │   ├── P512REPS20RES10240I50.txt
│   │   ├── P512REPS20RES1024I50.txt
│   │   └── P512REPS20RES23170I50.txt
│   └── SENDRECV
│       ├── P1024REPS20RES10240I50.txt
│       ├── P1024REPS20RES1024I50.txt
│       ├── P1024REPS20RES32768I50.txt
│       ├── P1REPS20RES10240I50.txt
│       ├── P1REPS20RES1024I50.txt
│       ├── P256REPS20RES10240I50.txt
│       ├── P256REPS20RES1024I50.txt
│       ├── P256REPS20RES16384I50.txt
│       ├── P32REPS20RES10240I50.txt
│       ├── P32REPS20RES1024I50.txt
│       ├── P32REPS20RES5793I50.txt
│       ├── P512REPS20RES10240I50.txt
│       ├── P512REPS20RES1024I50.txt
│       └── P512REPS20RES23170I50.txt
├── figures
│   ├── ALLTOALL
│   │   ├── StrongScalingExperiment
│   │   │   ├── PARALLELEFFICIENCY_REPS20I50.svg
│   │   │   ├── RUNTIME_REPS20RES10240I50.svg
│   │   │   ├── RUNTIME_REPS20RES1024I50.svg
│   │   │   └── SPEEDUP_REPS20I50.svg
│   │   └── WeakScalingExperiment
│   │       └── RUNTIME_REPS20RES1024I50.svg
│   └── SENDRECV
│       ├── StrongScalingExperiment
│       │   ├── PARALLELEFFICIENCY_REPS20I50.svg
│       │   ├── RUNTIME_REPS20RES10240I50.svg
│       │   ├── RUNTIME_REPS20RES1024I50.svg
│       │   └── SPEEDUP_REPS20I50.svg
│       └── WeakScalingExperiment
│           └── RUNTIME_REPS20RES1024I50.svg
└── src
    ├── PlotStrongScaling.py
    ├── PlotWeakScaling.py
    ├── arguments.hpp
    ├── game_sequential.hpp
    ├── main.cpp
    ├── mainATA.cpp
    ├── solver.hpp
    └── solverATA.hpp