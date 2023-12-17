# requirements on ubuntu
# sudo apt-get build-essentials
# sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev

# required modules on cluster
# module load mpi/openmpi-x86_64
# module load pmi/pmix-x86_64

CXX=/usr/local/bin/g++-13
MPICXX?=mpic++
CXXFLAGS := $(CXXFLAGS) -std=c++14 -O3 -Wall -pedantic -march=native

NPROCS?=$$((1)) $$((32)) $$((8*32)) $$((16*32)) $$((32*32))
RESOLUTION?=$$((1024*1024)) $$((10240*10240))

ITER?=50
REPETITION?=10

StrongScalingExperiment: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	@for RES in $(RESOLUTION) ; do \
		for NPROC in $(NPROCS) ; do \
			mpirun -n $$((NPROC)) --use-hwthread-cpus ./bin/GameOfLifeMPI $$((NPROC)) $(REPETITION) $$((RES)) $(ITER) ; \
		done ; \
	done

WeakScalingExperiment: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	@for RES in $$((1024*1024)) ; do \
		for NPROC in $(NPROCS) ; do \
			mpirun -n $$((NPROC)) --use-hwthread-cpus ./bin/GameOfLifeMPI $$((NPROC)) $(REPETITION) $$((RES*NPROC)) $(ITER) ; \
		done ; \
	done

DebugRun: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	mpirun -n 4 --use-hwthread-cpus ./bin/GameOfLifeMPI 4 $(REPETITION) 10 $(ITER)

clean:
	rm GameOfLifeMPI

