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
RESOLUTION?=$(1024) $(10240)

NPROCS_SMALL?=1 2 4
RESOLUTION_SMALL?=500 1000 1500

ITER?=20
REPETITION?=10

StrongScalingExperiment: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	@for RES in $(RESOLUTION) ; do \
		RES_ := $$((RES)) ; \
		for NPROC in $(NPROCS) ; do \
			NPROC_ = $$((NPROC)) ; \
			mpirun -n $(NPROC_) --use-hwthread-cpus ./bin/GameOfLifeMPI $(REPETITION) $(RES_) $(ITER) ; \
		done ; \
	done

StrongScalingExperimentSmall: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp ./src/PlotStrongScaling.py
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	@for RES in $(RESOLUTION_SMALL) ; do \
		for NPROC in $(NPROCS_SMALL) ; do \
			echo "Strong Scaling Experiment with Resolution of $$RES and $$NPROC Processors" ; \
			# mpirun -n $$NPROC --use-hwthread-cpus ./bin/GameOfLifeMPI $(REPETITION) $$RES $(ITER) ; \
		done ; \
	done
	python3 ./src/PlotStrongScaling.py --P $(NPROCS_SMALL) --REPS $(REPETITION) --RES $(RESOLUTION_SMALL) --I $(ITER)

WeakScalingExperiment: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	@for RES in $$((1024*1024)) ; do \
		for NPROC in $(NPROCS) ; do \
			mpirun -n $$((NPROC)) --use-hwthread-cpus ./bin/GameOfLifeMPI $(REPETITION) $$((RES*NPROC)) $(ITER) ; \
		done ; \
	done

WeakScalingExperimentSmall: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	@for RES in 100 ; do \
		for NPROC in $(NPROCS_SMALL) ; do \
			echo "Weak Scaling Experiment with Resolution of $$RES and $$NPROC Processors \n" ; \
			mpirun -n $$NPROC --use-hwthread-cpus ./bin/GameOfLifeMPI $(REPETITION) $$((RES*NPROC)) $(ITER) ; \
		done ; \
	done
	python3 ./src/PlotWeakScaling.py --P $(NPROCS_SMALL) --REPS $(REPETITION) --RES 100 --I $(ITER)

DebugRun: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	mpirun -n 4 --use-hwthread-cpus ./bin/GameOfLifeMPI $(REPETITION) 10 $(ITER)

clean:
	rm GameOfLifeMPI

