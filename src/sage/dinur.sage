# from mob import mob_transform
from mob_new import mob_transform
from fes import bruteforce
from fes import hamming_weight
import numpy as np
from random import randint
from utils import index_of, convert
from math import ceil
from fes_rec import fes_recover
import time

from collections import defaultdict

# For testing purposes
# _m_low = 5
# _m_high = 10
# _n_low = 2
# _n_high = 10

_m_low = 15
_m_high = 15 
_n_low = 15
_n_high = 15 

_time_bruteforce = 0
_time_u_values = 0
_time_mobius = 0
_time_full_eval = 0
_time_fetch_sol = 0
_time_output_potentials = 0
_time_solve_trials = 0
_time_solve = 0
_time_fes_recovery = 0

def random_system():
    system = []
    m = randint(_m_low, _m_high)
    n = randint(_n_low, _n_high)
    R = GF(2)[", ".join(["x" + str(i) for i in range(n)])]
    rem = 0
    for _ in range(m):
        f = R(GF(2)[R.gens()].random_element(degree=2, terms=Infinity))
        if (f != GF(2)(0)) and (f != GF(2)(1)):
            system.append(f)
        else:
            rem += 1
    m -= rem
    n1 = randint(1, n - 1)
    return system, n, m, n1, R

def random_system_with_sol():
    system, n, _, _, R = random_system()
    sol = randint(0, 2^n - 1) 
    system = [f if f(*convert(sol, n)) == 0 else (f + 1) for f in system]
    print("Forced solution: [" + ", ".join((f"{sol:0{n}b}")[::-1]) + "]")
    return system, sol, R, n

def test_u_values(trials, verbose=False):
    for _ in range(trials):
        system, n, m, n1, ring = random_system() 
        F_tilde = product((GF(2)(1) + f) for f in system)
        w = F_tilde.degree() - n1
        if verbose:
            print(f"Checking system\n{system}\n{(n, m, n1, w)}")
        V, ZV = compute_u_values(system, ring, n1, w, False)
        for y in range(2^(n - n1)):
            s0 = sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1)) for z_hat in range(2^n1))
            if V[y] != s0:
                print(f"Error found in V[{y}]\n\t{V[y]}\n\t{s0}")
                return y
            for i in range(n1):
                si = sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1 - 1)[:i], 0, *convert(z_hat, n1 - 1)[i:]) for z_hat in range(2^(n1 - 1)))
                if ((ZV[y] >> i) & 1) != si:
                    print(f"Error found in ZV[{i}][{y}]\n\t{((ZV[y] >> i) & 1)}\n\t{si}")
                    return y
        if verbose: 
            print("No errors found for system")
    print(f"No errors found in {trials} trial(s)")

def test_output_sol(trials, verbose=False):
    for _ in range(trials):
        system, n, m, n1, ring = random_system()
        F_tilde = product((GF(2)(1) + f) for f in system)
        w = F_tilde.degree() - n1
        if verbose:
            print(f"Checking system\n{system}\n{(n, m, n1, w)}")
        potentials = output_potentials(system, ring, n1, w, False)
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

# Not fully updated yet
def test_sliced_recovery(trials):
    for _ in range(trials):
        system, n, _, n1, R = random_system() 
        F_tilde = product((GF(2)(1) + f) for f in system)
        w = F_tilde.degree() - n1
        n = len(R.gens())
        R_sub = GF(2)[", ".join([str(var) for var in R.gens()[:n - n1]])]
         
        V, ZV, ZV_ = compute_u_values(system, R, n1, w + 1)
        print(ZV)
        print(ZV_)
        print("Computed U values")
        
        print("FES-based recovery")
        evals_0 = fes_recover(ZV, n - n1, w)
        print("Evaluation:") 
        print([stringify_bits(e, n1) for e in evals_0])
        
        print("Mob transform")
        U = np.full(n1 + 1, GF(2)(0))
        evals_1 = np.full((n1 + 1, 2^(n - n1)), GF(2)(0))
        # U[0] = mob_transform(V, R_sub.gens(), w)
        for i in range(1, n1 + 1):
            U[i] = mob_transform(ZV_[i - 1], R_sub.gens(), w)
        print("U_i", U)
        for i in range(n1 + 1):
            for y in range(2^(n - n1)):
                if U[i] in [0, 1]:
                    evals_1[i][y] = U[i]
                else:
                    evals_1[i][y] = U[i](*convert(y, n - n1))
        print(evals_1[1:])
        print("Finished mob and eval")

def dry_run_solve(trials):
    global _time_solve
    for i in range(trials):
        print(f"\n== {i} ==")
        system, known_sol, R, n = random_system_with_sol()
        print("Verify known solution:", [f(*convert(known_sol, n)) for f in system]) # Yes, I trust myself this "much"
        print("Known solution:", known_sol)
        print("System:", system, R)
        _time_solve -= time.time() 
        sol = solve(system, R)
        _time_solve += time.time()
        print("Solution found:", sol)
        print("Verify found solution:", [f(*sol) for f in system])

def GrayToBinary(num):
    mask = num
    while mask:
        mask >>= 1
        num ^^= mask
    return num

def compute_u_values(system, R, n1, w, slice):
    global _time_bruteforce
    global _time_u_values

    n = len(R.gens())

    _time_bruteforce -= time.time() 
    sols = bruteforce(system, R, n1, w + 1, slice)
    _time_bruteforce += time.time()

    V, ZV = None, None

    if slice:
        ZV = defaultdict(lambda: None)

        for i in range(2^(n - n1)):
          if hamming_weight(i) > w+1:
            #print(bin(i))
            continue
          ZV[i ^^ (i >> 1)] = 0
    else:
        V = defaultdict(lambda: GF(2)(0))
        ZV = [defaultdict(lambda: GF(2)(0)) for _ in range(n1)]

    _time_u_values -= time.time()

    for s in sols:
        y, z = s[:n - n1], s[n - n1:]

        idx = index_of(y)

        if slice:
          gidx = GrayToBinary(idx)
          yw = w + 1  # since we also interpolate U[0] for w+1 during fes_recover (bitsliced),
                      # we also need the (w+1)-weight interpolation point.
        else:
          gidx = idx
          yw = w

        if hamming_weight(gidx) <= yw:
            # print(y)
            idx = index_of(y)

            if slice:
                ZV[idx] ^^= 1
            else:
                V[idx] += GF(2)(1)

        if hamming_weight(gidx) <= w + 1:
          for i in range(1, n1 + 1):
            if z[i - 1] == 0:
                idx = index_of(y)

                if slice:
                    ZV[idx] ^^= (1 << i)
                else:
                    ZV[i - 1][idx] += GF(2)(1)

    _time_u_values += time.time()
    
    if slice:
        #print("ZV: ", ZV.keys())
        return ZV
    else:
        return V, ZV

def output_potentials(system, R, n1, w, fes_recovery):
    global _time_mobius
    global _time_full_eval
    global _time_fetch_sol
    global _time_fes_recovery

    n = len(R.gens())
    R_sub = GF(2)[", ".join([str(var) for var in R.gens()[:n - n1]])]

    #if fes_recovery:
    #  #VZV = compute_u_values(system, R, n1, w, fes_recovery)
    #  VZV = compute_u_values(system, R, n1, w+2, fes_recovery)
    #  #VZV = compute_u_values(system, R, n1, n-n1-1, fes_recovery)
    #  ## allVZV = compute_u_values(system, R, n1, n-n1, fes_recovery)

    #  ## for i in VZV.keys():
    #  ##   if VZV[i] != allVZV[i]:
    #  ##     print("err", i, bin(VZV[i]), bin(allVZV[i]))
    #  ##     VZV[i] ^^= 1
    #  ## #    VZV[i] = allVZV[i]


    #else:
    #  VZV = compute_u_values(system, R, n1, w, fes_recovery)

    VZV = compute_u_values(system, R, n1, w, fes_recovery)
     
    out = defaultdict(lambda: GF(2)(0))

    if fes_recovery:
        
        _time_fes_recovery -= time.time()
        
        evals = fes_recover(VZV, n - n1, w + 1)

        _time_fes_recovery += time.time()

        _time_fetch_sol -= time.time()

        for y_hat in range(2^(n - n1)):
            if evals[y_hat] & 1 == 1:

                if y_hat not in out: out[y_hat] = np.full(n1 + 1, GF(2)(0))

                out[y_hat][0] = GF(2)(1)
                for i in range(1, n1 + 1):
                    out[y_hat][i] = GF(2)((evals[y_hat] >> i) & 1) + 1# Bitsliced indexing

        _time_fetch_sol += time.time()
    else:
        U = []

        V, ZV = VZV

        _time_mobius -= time.time() 

        U.append(mob_transform(V, R_sub.gens(), w))
        for i in range(1, n1 + 1):
            U.append(mob_transform(ZV[i - 1], R_sub.gens(), w+1))

        _time_mobius += time.time()

        evals = np.full((n1 + 1, 2^(n - n1)), GF(2)(0))

        _time_full_eval -= time.time()


        evals = [None] * (n1 + 1)

        for i in range(n1 + 1):
            tmp = [0] * 2^(n-n1)
            for m in U[i].monomials():
              if m == 1:
                v = 0
              else:
                v = str(m)
                v = v.replace('x', '')
                v = [int(i) for i in v.split("*")]
                v = sum([2^i for i in v])
              tmp[v] = GF(2)(1)
    
            tmp = mob_transform(tmp, R_sub.gens())
    
            evals[i] = [0] * 2^(n-n1)
            for m in tmp.monomials():
              if m == 1:
                v = 0
              else:
                v = str(m)
                v = v.replace('x', '')
                v = [int(i) for i in v.split("*")]
                v = sum([2^i for i in v])
              evals[i][v] = GF(2)(1)

        ##print(evals[1], ZV[0])
        #for i in range(2^(n-n1)):
        #  if evals[0][i] != V[i]:
        #    print("err,", i, evals[1][i], ZV[0][i])

        #  for j in range(n1):
        #    if evals[j+1][i] != ZV[j][i]:
        #      print("err,", i, evals[j+1][i], ZV[j][i])


#        for i in range(n1 + 1):
#            for y in range(2^(n - n1)):
#                evals[i][y] = U[i](*convert(y, n - n1))

        _time_full_eval += time.time()

        out = defaultdict(lambda: GF(2)(0))

        _time_fetch_sol -= time.time()

        for y_hat in range(2^(n - n1)):
            if evals[0][y_hat] == 1:

                if y_hat not in out: out[y_hat] = np.full(n1 + 1, GF(2)(0))

                out[y_hat][0] = GF(2)(1)
                for i in range(1, n1 + 1):
                    out[y_hat][i] = evals[i][y_hat] + 1

        _time_fetch_sol += time.time()

    print("num pot sol: ", len(out))
    return out

def test_solution(system, sol):
    return not any(f(*sol) for f in system)

def preprocess(system, R):
    new_sys = []
    for f in system:
        for v in R.gens():
            if v^2 in f:
                f = f + v^2 + v
        new_sys.append(f)
    return new_sys

def gen_matrix_rank_l(l, m):
    A = random_matrix(GF(2), l, m)
    while A.rank() != l:
        A = random_matrix(GF(2), l, m)
    return A

def solve(system, R):
    global _time_solve_trials
    global _time_output_potentials
    system = preprocess(system, R)

    n = len(R.gens())
    n1 = int(ceil(n/(5.4))) # Quadratic systems are assumed here, see page 19 of full dinur paper for explanation
    l = n1 + 1
    m = len(system)

    sols = bruteforce(system, R, n, n)
    print("All solutions:", sols)

    potentials_solutions = []
    k = 0

    while k < 16:
        print("Commencing round", k)
        A = np.rint(np.random.rand(l, m))

        E_k = [sum(GF(2)(A[i][j]) * system[j] for j in range(m)) for i in range(l)]
        w = sum(f.degree() for f in E_k) - n1 # Do this differently 

        print("m: ", m, len(system))
        print("n: ", n)
        print("l: ", l, len(E_k))
        print("w: ", w)
        print("n1: ", n1)

        _time_output_potentials -= time.time()

        curr_potential_sol = output_potentials(E_k, R, n1, w, True)
        curr_potential_sol_mob = output_potentials(E_k, R, n1, w, False)
        #curr_potential_sol = output_potentials(E_k, R, n1, w, False)

        if len([f"{v:0{n}b}"[::-1] for v in curr_potential_sol.keys() if v not in curr_potential_sol_mob.keys()]) == 0 and \
           len([f"{v:0{n}b}"[::-1] for v in curr_potential_sol_mob.keys() if v not in curr_potential_sol.keys()]) == 0:
          print("OK - FES and Möbius agree!")

        _time_output_potentials += time.time()

        potentials_solutions.append(curr_potential_sol)

        _time_solve_trials -= time.time()
        for y_hat, potential_sol in curr_potential_sol.items(): # Iterate through solutions instead of all possible inputs
            for k1 in range(k):
                if all(potential_sol == potentials_solutions[k1][y_hat]):
                    sol = convert(y_hat, n - n1) + list(potential_sol[1:])
                    if test_solution(system, sol):
                        _time_solve_trials += time.time()
                        return sol
                    break
        _time_solve_trials += time.time()
        k += 1
    return None

def main():
    rounds = 1
    # test_u_values(1, False)
    # test_output_sol(1, True)
    # test_sliced_recovery(10)
    dry_run_solve(rounds)
    print("Bruteforce time:", _time_bruteforce/rounds)
    print("U value computation time:", _time_u_values/rounds)
    print("Time for FES interpolation:", _time_fes_recovery/rounds)
    print("Möbius transform time", _time_mobius/rounds)
    print("Full evaluation time", _time_full_eval/rounds)
    print("Fetch solutions time", _time_fetch_sol/rounds)
    print("Output potentials time", _time_output_potentials/rounds)
    print("Solve trials time:", _time_solve_trials/rounds)
    print("Solve time", _time_solve/rounds)

if __name__ == "__main__":
    main()
