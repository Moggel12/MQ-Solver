from mob import mob_transform
from fes import bruteforce
import numpy as np
from random import randint
from utils import index_of, convert
from math import ceil

# For testing purposes
_m_low = 5
_m_high = 10
_n_low = 2
_n_high = 10

def random_system():
    system = []
    m = randint(_m_low, _m_high)
    n = randint(_n_low, _n_high)
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
    return system, n, m, n1, R

def random_system_with_sol():
    system, n, _, _, R = random_system()
    sol = 2^n - 2 
    system = [f if f(*convert(sol, n)) == 0 else (f + 1) for f in system]
    return system, sol, R, n

def test_u_values(trials, verbose=False):
    for _ in range(trials):
        system, n, m, n1, ring = random_system() 
        F_tilde = product((GF(2)(1) + f) for f in system)
        w = F_tilde.degree() - n1
        if verbose:
            print(f"Checking system\n{system}\n{(n, m, n1, w)}")
        V, ZV = compute_u_values(system, R, n1, w)
        for y in range(2^(n - n1)):
            s0 = sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1)) for z_hat in range(2^n1))
            if V[y] != s0:
                print(f"Error found in V[{y}]\n\t{V[y]}\n\t{s0}")
                return y
            for i in range(n1):
                si = sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1 - 1)[:i], 0, *convert(z_hat, n1 - 1)[i:]) for z_hat in range(2^(n1 - 1)))
                if ZV[i][y] != si:
                    print(f"Error found in ZV[{i}][{y}]\n\t{ZV[i][y]}\n\t{si}")
                    return y
        if verbose: 
            print("No errors found for system")
    print(f"No errors found in {trials} trial(s)")

def test_output_sol(trials, verbose=False):
    for _ in range(trials):
        system, n, m, n1, R = random_system()
        F_tilde = product((GF(2)(1) + f) for f in system)
        w = F_tilde.degree() - n1
        if verbose:
            print(f"Checking system\n{system}\n{(n, m, n1, w)}")
        potentials = output_potentials(system, R, n1, w)
        for y in range(2^(n - n1)):
            s0 = sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1)) for z_hat in range(2^n1))
            if s0 == GF(2)(1):
                if y not in potentials:
                    print(f"Error in out: {y} is not a key")
                    return y
                if potentials[y][0] != GF(2)(1):
                    print(f"Error in out[{y}][0]: Not equal to 1 when sum evaluates to 1")
                    return y
                for i in range(1, n1 + 1):
                    si = sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1 - 1)[:i - 1], 0, *convert(z_hat, n1 - 1)[i - 1:]) for z_hat in range(2^(n1 - 1)))
                    if potentials[y][i] == si:
                        print(f"Error in out[{y}][{i}]: Value NOT inverse of the evaluation of U_{i}({y}) = {si} = {potentials[y][i]}")
                        return y
        if verbose:
            print("No errors found for system")
    print(f"No errors found in {trials} trial(s)")

def dry_run_solve(trials):
    for _ in range(trials):
        print("\n===")
        system, known_sol, R, n = random_system_with_sol()
        print("Verify known solution:", [f(*convert(known_sol, n)) for f in system]) # Yes, I trust myself this "much"
        print("Known solution:", known_sol)
        print("System:", system, R)
        sol = solve(system, R)
        print("Solution found:", sol)
        print("Verify found solution:", [f(*sol) for f in system])

def compute_u_values(system, R, n1, w):
    n = len(R.gens())
    sols = bruteforce(system, R, n1, w + 1)
    l = [math.comb(n - n1, i) for i in range(w + 2)]
    V = {i: GF(2)(0) for i in range(sum(l[:-1]))}
    ZV = {i: {j: GF(2)(0) for j in range(sum(l))} for i in range(n1)}
    for s in sols:
        y, z = s[:n - n1], s[n - n1:]
        if sum(y) <= w:
            idx = index_of(y)
            V[idx] += 1
        for i in range(n1):
            if z[i] == 0:
                idx = index_of(y)
                ZV[i][idx] += 1
    return V, ZV

def output_potentials(system, R, n1, w):
    n = len(R.gens())
    R_sub = GF(2)[", ".join([str(var) for var in R.gens()[:n - n1]])]
    V, ZV = compute_u_values(system, R, n1, w + 1)
    U = np.full(n1 + 1, GF(2)(0))
    U[0] = mob_transform(V, R_sub.gens())
    for i in range(1, n1 + 1):
        U[i] = mob_transform(ZV[i - 1], R_sub.gens())
    evals = np.full((n1 + 1, 2^(n - n1)), GF(2)(0))
    for i in range(n1 + 1):
        for y in range(2^(n - n1)):
            evals[i][y] = U[i](*convert(y, n - n1))
    out = np.full((2^(n - n1), n1 + 1), GF(2)(0))
    for y_hat in range(2^(n - n1)):
        if evals[0][y_hat] == 1:
            out[y_hat][0] = GF(2)(1)
            for i in range(1, n1 + 1):
                out[y_hat][i] = evals[i][y_hat] + 1
    return out

def test_solution(system, sol):
    return not any(f(*sol) for f in system)

def solve(system, R):
    n = len(R.gens())
    n1 = int(ceil(n/(5.4))) # Quadratic systems are assumed here, see page 19 of full dinur paper for explanation
    l = n1 + 1
    m = len(system)
    potentials_solutions = []
    k = 0
    while True:
        print("Commencing round", k)
        A = np.rint(np.random.rand(l, m))
        E_k = [sum(GF(2)(A[i][j]) * system[j] for j in range(m)) for i in range(l)]
        w = product((1 + f) for f in system).degree() - n1 # Do this differently
        curr_potential_sol = output_potentials(system, R, n1, w)
        potentials_solutions.append(curr_potential_sol)
        for y_hat in range(2^(n - n1)):
            if curr_potential_sol[y_hat][0] == 1:
                for k1 in range(k):
                    if all(curr_potential_sol[y_hat] == potentials_solutions[k1][y_hat]):
                        sol = convert(y_hat, n - n1) + list(curr_potential_sol[y_hat][1:])
                        if test_solution(system, sol):
                            return sol
                        break
        k += 1

def main():
    # test_u_values(1, False)
    # test_output_sol(20, False)
    dry_run_solve(10)

if __name__ == "__main__":
    main()
