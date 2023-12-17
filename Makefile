# requirements on ubuntu
# sudo apt-get build-essentials
# sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev

# required modules on cluster
# module load mpi/openmpi-x86_64
# module load pmi/pmix-x86_64

CXX=/usr/local/bin/g++-13
MPICXX?=mpic++
CXXFLAGS := $(CXXFLAGS) -std=c++14 -O3 -Wall -pedantic -march=native -ffast-math
NPROCS=4
RESOLUTION=10
ITER=5
REPETITION=10


GameOfLifeMPI: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)

run: Makefile ./bin/GameOfLifeMPI
	mpirun -n $(NPROCS) --use-hwthread-cpus ./bin/GameOfLifeMPI $(NPROCS) $(REPETITION) $(RESOLUTION) $(ITER)

clean:
	rm GameOfLifeMPI

