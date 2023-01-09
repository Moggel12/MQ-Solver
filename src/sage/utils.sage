import numpy as np
from itertools import combinations
import math

def bitslice(f_sys, vars):
    f_sys_sliced = np.zeros(math.comb(len(vars), 2) + len(vars) + 1, dtype=int)
    for j, poly in enumerate(f_sys):
        if poly in [GF(2)(0), GF(2)(1)]:
            f_sys_sliced[0] ^^= int(poly) << j
        else:
            f_sys_sliced[0] ^^= int(poly.constant_coefficient()) << j
            i = 1
            for v in vars:
                f_sys_sliced[i] ^^= int(poly.coefficient({v_: 1 if v == v_ else 0 for v_ in vars})) << j
                i += 1
            for v1, v2 in combinations(range(len(vars)), 2):
                f_sys_sliced[i] ^^= int(poly.coefficient({vars[v1]: 1, vars[v2]: 1, **{v: 0 for v in vars if v not in [vars[v1], vars[v2]]}})) << j
                i += 1
    return f_sys_sliced

def convert(v, n):
    v = bin(v)[2:]
    return list(map(lambda i : (int(v[-i]) if i <= len(v) else 0), range(1, n + 1)))

def index_of(y_list):
    return sum(b << i for i, b in enumerate(y_list))
