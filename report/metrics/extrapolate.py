#! /usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit
import os
import re
import sys
import pandas as pd

def _exp(x, m, t, b):
    return m * np.exp2(t * x) + b

def check_file(fname, start, stop):
    return (
            not fname.startswith("procinfo")
        ) or (
            start > int(fname.split("_")[-1]) 
        ) or (
            int(fname.split("_")[-1]) > stop
        )


def mem():
    REGEX = r"VmHWM:\s+(\d+) kB"

    data = {}

    if len(sys.argv) != 4:
        print("No arguments given")
        sys.exit(1)

    directory = sys.argv[2]

    range_str = sys.argv[3]

    range_val = [int(v) for v in range_str.split("..")]

    for fname in os.listdir(directory):
        if check_file(fname, range_val[0], range_val[1]): continue
        n = int(fname.split("_")[-1])
        with open(os.path.join(directory, fname), "r") as f:
            m = re.search(REGEX, f.read())
            if m:
                data[n] = int(m.group(1))

    xs = np.array(sorted(list(data.keys())), dtype=np.float64)
    ys = np.array([data[x] for x in xs], dtype=np.float64)

    p0 = (1, 1, 0)
    params, cv = curve_fit(_exp, xs, ys, p0, bounds=(0, [np.inf, np.inf, np.inf]))
    m, t, b = params

    squared_diffs = np.square(ys - _exp(xs, m, t, b))
    squared_diffs_from_mean = np.square(ys - np.mean(ys))
    r_squared = 1 - np.sum(squared_diffs) / np.sum(squared_diffs_from_mean)

    print(f"R^2 = {r_squared}")
    print(f"\\addplot [] {{{m} * 2^({t} * x) + {b}}};")
    print(f"\\addplot [] coordinates {{{''.join(str(t) for t in zip(xs, ys))}}};")


def timings():
    if len(sys.argv) != 3:
        print("No arguments given")
        sys.exit(1)

    fname = sys.argv[2]

    data = pd.read_csv(fname)

    xs = data["n"]

    y_std = np.array(data["standard"], dtype=np.float64)
    y_vec = np.array(data["vectorized"], dtype=np.float64)
    y_fes = np.array(data["fes"], dtype=np.float64)

    lim = (0,30)
    x_lim = np.linspace(*lim,100)

    for i, col in enumerate(data):
        if col == "n": continue
        print(f"#### {col} ####")
        ys = np.array(data[col], dtype=np.float64)

        p0 = (1, 1, 0)
        params, cv = curve_fit(_exp, xs, ys, p0,bounds=(0, [np.inf, np.inf, np.inf]))
        m, t, b = params

        squared_diffs = np.square(ys - _exp(xs, m, t, b))
        squared_diffs_from_mean = np.square(ys - np.mean(ys))
        r_squared = 1 - np.sum(squared_diffs) / np.sum(squared_diffs_from_mean)
        print(f"R^2 = {r_squared}")
        
        print(f"\\addplot [] {{{m} * 2^({t} * x) + {b}}};")

if __name__ == "__main__":
    if (len(sys.argv) < 3) or (len(sys.argv) > 4):
        print("Not enough arguments")
        sys.exit(1)
    ext_from = sys.argv[1]
    if ext_from == "mem":
        mem()
    elif ext_from == "time":
        timings()
    else:
        print("Invalid choice of first arugment")
        sys.exit(1)