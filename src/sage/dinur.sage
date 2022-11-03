from mob import mob_transform
from fes import bruteforce, convert
import numpy as np
from random import randint
# from monotonic_gray import monotonic_bounded

def index_of(y_list):
    return sum(b * 2^i for i, b in enumerate(y_list))

def test_u_values(trials, verbose=False):
    for _ in range(trials):
        system = []
        m = randint(5, 10)
        n = randint(2, 10)
        R = GF(2)[", ".join(["x" + str(i) for i in range(n)])]
        rem = 0
        for _ in range(m):
            f = R(GF(2)[R.gens()].random_element(degree=2))
            if (f != GF(2)(0)) and (f != GF(2)(1)):
                system.append(f)
            else:
                rem += 1
        m -= rem
        n1 = randint(1, n - 1)
        F_tilde = product((GF(2)(1) + f) for f in system)
        w = F_tilde.degree() - n1
        if verbose:
            print(f"Checking system\n{system}\n{(n, m, n1, w)}")
        V, ZV = compute_u_values(system, R, n1, w)
        for y in range(2^(n - n1)):
            if V[y] != sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1)) for z_hat in range(2^n1)):
                print(f"Error found in V[{y}]\n\t{V[y]}\n\t{sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1)) for z_hat in range(2^n1))}")
                return y
            for i in range(n1):
                if ZV[i][y] != sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1 - 1)[:i], 0, *convert(z_hat, n1 - 1)[i:]) for z_hat in range(2^(n1 - 1))):
                    print(f"Error found in ZV[{i}][{y}]\n\t{ZV[i][y]}\n\t{sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1 - 1)[:i], 0, *convert(z_hat, n1 - 1)[i:]) for z_hat in range(2^(n1 - 1)))}")
                    return y
        if verbose: 
            print("No errors found for system")
    print(f"No errors found in {trials} trials")

def test_output_sol():
    pass

def compute_u_values(system, R, n1, w):
    n = len(R.gens())
    sols = bruteforce(system, R, n1, w + 1)
    l = [math.comb(n - n1, i) for i in range(w + 2)]
    V = {i: GF(2)(0) for i in range(sum(l[:-1]))}
    ZV = {i: {j: GF(2)(0) for j in range(sum(l))} for i in range(n1)}
    for s in sols:
        y,z = s[:n - n1], s[n - n1:]
        if sum(y) <= w:
            idx = index_of(y)
            V[idx] += 1
        for i in range(n1):
            if z[i] == 0:
                idx = index_of(y)
                ZV[i][idx] += 1
    return V, ZV

def output_potentials(system, R, n1, w):
    V, ZV = compute_u_values(system, n, n1, w + 1)
    U = np.full(n1 + 1, GF(2)(0))
    U[0] = mob_transform(V, R.gens()[:n - n1])
    for i in range(1, n1 + 1):
        U[i] = mob_transform(ZV[i - 1], R.gens()[:n - n1])
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

def main():
    pass

if __name__ == "__main__":
    test_u_values(100)