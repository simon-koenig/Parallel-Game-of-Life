# requirements on ubuntu
# sudo apt-get build-essentials
# sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev

# required modules on cluster
# module load mpi/openmpi-x86_64
# module load pmi/pmix-x86_64

CXX=/usr/local/bin/g++-13
MPICXX?=mpic++
CXXFLAGS := $(CXXFLAGS) -std=c++14 -O3 -Wall -pedantic -march=native

NPROCS ?= 32 256 512 1024
RESOLUTION?=1024 10240

NPROCS_SMALL?=1 2 4
RESOLUTION_SMALL?=500 1000 1500

ITER?=10
REPETITION?=20

StrongScalingExperiment: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	@for RES in $(RESOLUTION) ; do \
		echo "Strong Scaling Experiment with Resolution of $$RES and 1 Processor"; \
		srun -t 5 -p q_student -N 1 --ntasks-per-node=1 ./bin/GameOfLifeMPI $(REPETITION) $$RES $(ITER); \
		for NPROC in $(NPROCS) ; do \
			echo "Strong Scaling Experiment with Resolution of $$RES and $$NPROC Processors" ; \
			srun -t 5 -p q_student -N $$((NPROC/32)) --ntasks-per-node=32 ./bin/GameOfLifeMPI $(REPETITION) $$RES $(ITER) ; \
		done ; \
	done

EvaluateStrongScalingExperiment: ./src/PlotStrongScaling.py
	python3 ./src/PlotStrongScaling.py --P $(NPROCS) --REPS $(REPETITION) --RES $(RESOLUTION) --I $(ITER)

StrongScalingExperimentSmall: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp ./src/PlotStrongScaling.py
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	@for RES in $(RESOLUTION_SMALL) ; do \
		for NPROC in $(NPROCS_SMALL) ; do \
			echo "Strong Scaling Experiment with Resolution of $$RES and $$NPROC Processors" ; \
			mpirun -n $$NPROC --use-hwthread-cpus ./bin/GameOfLifeMPI $(REPETITION) $$RES $(ITER) ; \
		done ; \
	done
	python3 ./src/PlotStrongScaling.py --P $(NPROCS_SMALL) --REPS $(REPETITION) --RES $(RESOLUTION_SMALL) --I $(ITER)

WeakScalingExperiment: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	@for RES in 1024 ; do \
		echo "Weak Scaling Experiment with 1 Processor and Resolution of $$RES"; \
		srun -t 5 -p q_student -N 1 --ntasks-per-node=1 ./bin/GameOfLifeMPI $(REPETITION) $$RES $(ITER); \
		for NPROC in $(NPROCS) ; do \
			tmp=$$( (echo "scale=20;(sqrt($$NPROC)*$$RES)"|bc) ); \
			tmp=$$( (printf "%.f" $$tmp) ); \
			echo "Weak Scaling Experiment with $$NPROC Processors and Resolution of $$tmp" ; \
			srun -t 5 -p q_student -N $$((NPROC/32)) --ntasks-per-node=32 ./bin/GameOfLifeMPI $(REPETITION) $$((tmp)) $(ITER) ; \
		done ; \
	done

EvaluateWeakScalingExperiment: ./src/PlotWeakScaling.py
	python3 ./src/PlotWeakScaling.py --P $(NPROCS) --REPS $(REPETITION) --RES 1024 --I $(ITER)	

WeakScalingExperimentSmall: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	for RES in 500 ; do \
		for NPROC in $(NPROCS) ; do \
			tmp=`awk ' { print sqrt($$NPROC) }'`;\
			echo "Weak Scaling Experiment with Resolution of $$RES and $${tmp} Processors \n";\
			mpirun -n $$NPROC --use-hwthread-cpus ./bin/GameOfLifeMPI $(REPETITION) $$((RES*NPROC)) $(ITER);\
		done ; \
	done
	python3 ./src/PlotWeakScaling.py --P $(NPROCS_SMALL) --REPS $(REPETITION) --RES 500 --I $(ITER)

DebugRun: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	mpirun -n 4 --use-hwthread-cpus ./bin/GameOfLifeMPI $(REPETITION) 11 $(ITER)

SingleProcessor: Makefile ./src/main.cpp ./src/solver.hpp ./src/arguments.hpp
	$(MPICXX) ./src/main.cpp -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	srun -t 5 -p q_student -N 1 --tasks-per-node=1 ./bin/GameOfLifeMPI $(REPETITION) 200 $(ITER)

clean:
	rm GameOfLifeMPI

