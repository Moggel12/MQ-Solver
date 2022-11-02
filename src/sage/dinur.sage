from mob import mob_transform
from fes import bruteforce
from math import comb
import numpy as np
from monotonic_gray import monotonic_bounded


def index_of(y_list):
    return sum([b * 2^i for i, b in enumerate(y_list)])

def compute_u_values(system, n, n1, w):
    g = subspace_gen(n, n1, w)
    sols = bruteforce(system, R, n - n1, w + 1, lambda _ : next(g)) # Not the prettiest code I've written
    l = [math.comb(n - n1, i) for i in range(w + 2)]
    V = {i: GF(2)(0) for i in range(sum(l[:-1]))}
    ZV = {i: {j: GF(2)(0) for j in range(sum(l))}  for i in range(n1)}
    for s in sols:
        y,z = sols[:n-n1], sols[n - n1:]
        if sum(y) <= w:
            idx = index_of(y)
            V[idx] += 1
        for i in range(1, n1 + 1):
            if z[i] == 0:
                idx = index_of(y)
                ZV[i - 1][idx] += 1 # NOTICE: zero-indexing => i - 1 (reminder for future bugs)
    return V, ZV

def output_potentials():
    pass

def solve():
    pass

def test_u_values(system, n1):
    

def main():
    pass

if __name__ == "__main__":
    pass