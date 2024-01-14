MPICXX?=mpic++
CXXFLAGS := $(CXXFLAGS) -std=c++14 -O3 -Wall -pedantic -march=native

NPROCS ?=32 256 512 1024
RESOLUTION?=1024 10240

ITER?=50
REPETITION?=20

TEST_CORRECTNESS?=0 # if set to 1, visual representations of parallel and sequential result will be printed
IMPLEMENTATION=?0 # Implementation 0 = SEND_RECV, Implementation 1 = ALLTOALL


SOLVER:=./src/solver.hpp
MAIN:=./src/main.cpp

ifeq ($(IMPLEMENTATION), 1)
SOLVER=./src/solverATA.hpp
MAIN=./src/mainATA.cpp
endif

StrongScalingExperiment: Makefile ./src/arguments.hpp $(MAIN) $(SOLVER)
	$(MPICXX) $(MAIN) -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	@for RES in $(RESOLUTION) ; do \
		echo "Strong Scaling Experiment with Resolution of $$RES and 1 Processor"; \
		srun -t 5 -p q_student -N 1 --ntasks-per-node=1 ./bin/GameOfLifeMPI $(REPETITION) $$RES $(ITER) $(TEST_CORRECTNESS); \
		for NPROC in $(NPROCS) ; do \
			echo "Strong Scaling Experiment with Resolution of $$RES and $$NPROC Processors" ; \
			srun -t 5 -p q_student -N $$((NPROC/32)) --ntasks-per-node=32 ./bin/GameOfLifeMPI $(REPETITION) $$RES $(ITER) $(TEST_CORRECTNESS); \
		done ; \
	done

EvaluateStrongScalingExperiment: ./src/PlotStrongScaling.py
	python3 ./src/PlotStrongScaling.py --P $(NPROCS) --REPS $(REPETITION) --RES $(RESOLUTION) --I $(ITER) --IMPL $(IMPLEMENTATION)

WeakScalingExperiment: Makefile $(MAIN) $(SOLVER) ./src/arguments.hpp
	$(MPICXX) $(MAIN) -o ./bin/GameOfLifeMPI -lpthread -DUSEMPI $(CXXFLAGS)
	@for RES in 1024 ; do \
		echo "Weak Scaling Experiment with 1 Processor and Resolution of $$RES"; \
		srun -t 5 -p q_student -N 1 --ntasks-per-node=1 ./bin/GameOfLifeMPI $(REPETITION) $$RES $(ITER) $(TEST_CORRECTNESS); \
		for NPROC in $(NPROCS) ; do \
			tmp=$$( (echo "scale=20;(sqrt($$NPROC)*$$RES)"|bc) ); \
			tmp=$$( (printf "%.f" $$tmp) ); \
			echo "Weak Scaling Experiment with $$NPROC Processors and Resolution of $$tmp" ; \
			srun -t 5 -p q_student -N $$((NPROC/32)) --ntasks-per-node=32 ./bin/GameOfLifeMPI $(REPETITION) $$((tmp)) $(ITER) $(TEST_CORRECTNESS); \
		done ; \
	done

EvaluateWeakScalingExperiment: ./src/PlotWeakScaling.py
	python3 ./src/PlotWeakScaling.py --P $(NPROCS) --REPS $(REPETITION) --RES 1024 --I $(ITER) --IMPL $(IMPLEMENTATION)	


CleanBinary:
	rm bin/GameOfLifeMPI

CleanData:
	rm data/ALLTOALL/*
	rm data/SENDRECV/*

CleanFigures:
	rm figures/ALLTOALL/StrongScalingExperiment/*
	rm figures/ALLTOALL/WeakScalingExperiment/*
	rm figures/SENDRECV/StrongScalingExperiment/*
	rm figures/SENDRECV/WeakScalingExperiment/*
	

