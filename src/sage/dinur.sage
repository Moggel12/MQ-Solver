from mob import mob_transform
from fes import bruteforce, convert
import numpy as np
# from monotonic_gray import monotonic_bounded

def index_of(y_list):
    return sum(b * 2^i for i, b in enumerate(y_list))

def compute_u_values(system, n, n1, w):
    sols = bruteforce(system, R, n - n1, w + 1)
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
                ZV[i - 1][idx] += 1
    return V, ZV

def output_potentials(system, R, n1, w):
    V, ZV = compute_u_values(system, n, n1, w + 1)
    U = np.full(n1 + 1, GF(2)(0))
    U[0] = mob_transform(V, R.gens()[:n - n1])
    for i in range(1, n1 + 1):
        U[i] = mob_transform(ZV[i], R.gens()[:n - n1])
    evals = np.full((n1 + 1, 2^(n - n1)), GF(2)(0))
    for i in range(n1 + 1):
        for y in range(2^(n - n1)):
            evals[i][y] = U[i](*convert(y, n1))
    out = np.full((2^(n - n1), n1 + 1), GF(2)(0))
    for y_hat in range(2^(n - n1)):
        if evals[0][y_hat] == 1:
            out[y_hat][0] = GF(2)(1)
            for i in range(1, n1 + 1):
                out[y][i] = evals[i][y] + 1
    return out

def solve():
    pass

def test_u_values(system, n1):
    pass

def main():
    pass

if __name__ == "__main__":
    pass