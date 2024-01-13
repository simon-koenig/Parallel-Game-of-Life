# Usual imports for plotting
import argparse
import os
from pathlib import Path
import itertools
import matplotlib.pyplot as plt
import numpy as np


def read_timing_file(filename):
    summed_time = []
    mean_time = []
    max_time = []
    with open(filename, "r") as file:
        for line in itertools.islice(file, 7, None):
            summed_time.append(float(line.split(",")[0]))
            mean_time.append(float(line.split(",")[1]))
            max_time.append(float(line.split(",")[2]))    
    return {"summed_time": summed_time, "mean_time": mean_time, "max_time": max_time}    


# Parsing command line arguments to get details about the benchmark.
parser = argparse.ArgumentParser(description='Supply number of processors, performed repetitions, used resolutions and performed iterations.')

parser.add_argument('--P', dest="P", type=int, nargs='+',
                    help='Supply number of processors benchmarked')
parser.add_argument('--REPS', dest="REPS", type=int, nargs='+',
                    help='Supply performed repetiotions of the benchmarks')
parser.add_argument('--RES', dest="RES", type=int, nargs='+',
                    help='Supply grid resolution of the benchmarks')
parser.add_argument('--I', dest="I", type=int, nargs='+',
                    help='Supply number of iterations of the benchmarks')
parser.add_argument('--IMPL', dest="IMPL", type=int,
                    help='Supply implementation which was used')
args = parser.parse_args()

# Check that exactly one number of repetitions and iterations is supplied.
assert len(args.REPS) == 1, "Plotting for multiple numbers of repetitions is not supported."
assert len(args.I) == 1, "Plotting for multiple numbers of iterations is not supported."

# Get data from right folder according to the used implementation
if (args.IMPL == 0):
    datapath = Path.cwd() / "data" / "SENDRECV"
    figurepath = Path.cwd() / "figures" / "SENDRECV"/ "StrongScalingExperiment"
elif (args.IMPL == 1):
    datapath = Path.cwd() / "data" / "ALLTOALL"
    figurepath = Path.cwd() / "figures" / "ALLTOALL" / "StrongScalingExperiment"

# Insert sequential run
args.P.insert(0, 1)

# Iterate over all relevant files and read in data.
data = {}
for combination in itertools.product(args.P, args.REPS, args.RES, args.I):
    filename = "/P" + str(combination[0]) + "REPS" + str(combination[1]) + "RES" + str(combination[2]) + "I" + str(combination[3]) + ".txt"
    data[combination] = read_timing_file(str(datapath) + filename)

speedup_plot, speedup_axis = plt.subplots()
pareff_plot, pareff_axis = plt.subplots()

for res in args.RES:
    runtime_plot, runtime_axis = plt.subplots()

    mean_rt = []
    std_rt = []
    mean_su = []
    std_su = []
    mean_pe = []
    std_pe = []

    # getting sequential runtime for calculating speedup and parallel efficiency
    rt_seq = np.mean(data[(1, args.REPS[0], res, args.I[0])]["max_time"])

    for p in args.P:

        mean_rt.append(np.mean(data[(p, args.REPS[0], res, args.I[0])]["max_time"]))
        std_rt.append(np.std(data[(p, args.REPS[0], res, args.I[0])]["max_time"]))

        if p != 1:

            mean_su.append(np.mean(rt_seq / data[(p, args.REPS[0], res, args.I[0])]["max_time"]))
            std_su.append(np.std(rt_seq / data[(p, args.REPS[0], res, args.I[0])]["max_time"]))
            mean_pe.append(np.mean(rt_seq / np.multiply(p, data[(p, args.REPS[0], res, args.I[0])]["max_time"])))
            std_pe.append(np.std(rt_seq / np.multiply(p, data[(p, args.REPS[0], res, args.I[0])]["max_time"])))

    # create and save individual runtime plot for each resolution
    runtime_axis.plot(args.P, mean_rt, "-")
    runtime_axis.set_yscale('log')
    runtime_axis.fill_between(args.P, np.subtract(mean_rt, std_rt), np.add(mean_rt, std_rt), alpha=0.2)
    runtime_axis.set_title("Mean runtimes of resolution " + str(res) + " for " + str(args.I[0]) + " iterations") 
    runtime_axis.set_xlabel("Number of processes")
    runtime_axis.set_ylabel("Mean runtime in sec")
    runtime_plot.savefig(str(figurepath) + "/RUNTIME_REPS" + str(args.REPS[0]) + "RES" + str(res) + "I" + str(args.I[0]) + ".svg")

    # add speedups of current resolution to plot
    speedup_axis.plot(args.P[1:], mean_su, "-", label="RES = " + str(res))
    speedup_axis.fill_between(args.P[1:], np.subtract(mean_su, std_su), np.add(mean_su, std_su), alpha=0.2)

    # add parallel efficiency of current resolution to plot
    pareff_axis.plot(args.P[1:], mean_pe, "-", label="RES = " + str(res))
    pareff_axis.fill_between(args.P[1:], np.subtract(mean_pe, std_pe), np.add(mean_pe, std_pe), alpha=0.2)

# Adapt and save speed-up plot
speedup_axis.set_title("Mean speed-up for " + str(args.I[0]) + " iterations") 
speedup_axis.set_xlabel("Number of processes")
speedup_axis.set_ylabel("Mean speed-up")
speedup_axis.legend()
speedup_plot.savefig(str(figurepath) + "/SPEEDUP_REPS" + str(args.REPS[0]) + "I" + str(args.I[0]) + ".svg")


# Adapt and save parallel efficiency plot
pareff_axis.set_title("Mean parallel efficiency for " + str(args.I[0]) + " iterations") 
pareff_axis.set_xlabel("Number of processes")
pareff_axis.set_ylabel("Mean parallel effiency")
pareff_axis.legend()
pareff_plot.savefig(str(figurepath) + "/PARALLELEFFICIENCY_REPS" + str(args.REPS[0]) + "I" + str(args.I[0]) + ".svg")
