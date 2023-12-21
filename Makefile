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

ITER?=20
REPETITION?=2

StrongScalingExperiment: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	@for RES in $(RESOLUTION) ; do \
		RES_ := $$((RES)) ; \
		for NPROC in $(NPROCS) ; do \
			NPROC_ = $$((NPROC)) ; \
			mpirun -n $(NPROC_) --use-hwthread-cpus ./bin/GameOfLifeMPI $(REPETITION) $(RES_) $(ITER) ; \
		done ; \
	done

WeakScalingExperiment: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	@for RES in $$((1024*1024)) ; do \
		for NPROC in $(NPROCS) ; do \
			mpirun -n $$((NPROC)) --use-hwthread-cpus ./bin/GameOfLifeMPI $(REPETITION) $$((RES*NPROC)) $(ITER) ; \
		done ; \
	done

DebugRun: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	mpirun -n 4 --use-hwthread-cpus ./bin/GameOfLifeMPI $(REPETITION) 10 $(ITER)

clean:
	rm GameOfLifeMPI

