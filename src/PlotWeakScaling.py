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
args = parser.parse_args()

# Check that exactly one number of repetitions and iterations is supplied.
assert len(args.REPS) == 1, "Plotting for multiple numbers of repetitions is not supported."
assert len(args.I) == 1, "Plotting for multiple numbers of iterations is not supported."

# Insert sequential run
args.P.insert(0, 1)

# Get path of "data"-folder and "figures"-folder
datapath = Path.cwd() / "data"
figurepath = Path.cwd() / "figures" / "WeakScalingExperiment"

# Iterate over all relevant files and read in data.
data = {}

for combination in itertools.product(args.P, args.REPS, args.RES, args.I):
    filename = "/P" + str(combination[0]) + "REPS" + str(combination[1]) + "RES" + str(round(combination[2]*np.sqrt(combination[0]))) + "I" + str(combination[3]) + ".txt"
    data[combination] = read_timing_file(str(datapath) + filename)

for res in args.RES:
    runtime_plot, runtime_axis = plt.subplots()

    mean_rt = []
    std_rt = []

    # getting sequential runtime for calculating speedup and parallel efficiency
    rt_seq = np.mean(data[(1, args.REPS[0], res, args.I[0])]["max_time"])

    for p in args.P:

        mean_rt.append(np.mean(data[(p, args.REPS[0], res, args.I[0])]["max_time"]))
        std_rt.append(np.std(data[(p, args.REPS[0], res, args.I[0])]["max_time"]))

    # create and save individual runtime plot for each resolution
    runtime_axis.plot(args.P, mean_rt, "-")
    runtime_axis.fill_between(args.P, np.subtract(mean_rt, std_rt), np.add(mean_rt, std_rt), alpha=0.2)
    runtime_axis.set_title("Mean runtimes of resolution " + str(res) + " (P=1) for " + str(args.I[0]) + " iterations") 
    runtime_axis.set_xlabel("Number of processes")
    runtime_axis.set_ylabel("Mean runtime in sec")
    runtime_plot.savefig(str(figurepath) + "/RUNTIME_REPS" + str(args.REPS[0]) + "RES" + str(res) + "I" + str(args.I[0]) + ".svg")









