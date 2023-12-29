# Usual imports for plotting
import argparse
import os
from pathlib import Path

# Parsing command line arguments to get details about the benchmark.
parser = argparse.ArgumentParser(description='Supply number of processors, performed repetitions, used resolutions and performed iterations.')

parser.add_argument('--P', type=int, nargs='+',
                    help='Supply number of processors benchmarked')
parser.add_argument('--REPS', type=int, nargs='+',
                    help='Supply performed repetiotions of the benchmarks')
parser.add_argument('--RES', type=int, nargs='+',
                    help='Supply grid resolution of the benchmarks')
parser.add_argument('--I', type=int, nargs='+',
                    help='Supply number of iterations of the benchmarks')
args = parser.parse_args()

datapath = Path.cwd() / "data"



