from mob_new import mob_transform
from fes import bruteforce
import numpy as np
from random import randint
from utils import index_of, convert, random_systems, fetch_systems_interactive, CLEAR, WARNING, FAIL, SUCCESS, ITER
from math import ceil
from fes_rec import fes_recover
import time

from c_config import Type, srand, rand, call_c_test, RSEED, C_POLY_T, C_VARS_T

import ctypes as ct
from utils import bitslice, fetch_c_func, write_fukuoka

from collections import defaultdict

_m_low  = 5 
_n_low  = 5
_m_high = 5 
_n_high = 5 

_time_bruteforce = 0
_time_u_values = 0
_time_mobius = 0
_time_full_eval = 0
_time_fetch_sol = 0
_time_output_potentials = 0
_time_solve_trials = 0
_time_solve = 0
_time_fes_recovery = 0

def test_c_output_solutions(sys_tuple):
    # system, _, ring, n = random_systems_with_sol(_m_low, _m_high, _n_low, _n_high)
    system, n, _, ring, _ = sys_tuple
    n1 = int(ceil(n/(5.4))) # Quadratic systems are assumed here, see page 19 of full dinur paper for explanation
    w = sum(f.degree() for f in system) - n1
    print(system)
    print(n, n1, w)
    sl_system = bitslice(system, ring.gens())
    c_solutions = c_output_potentials(sl_system, n, n1, w)
    if c_solutions == None:
        print(f"{FAIL}Error with memory allocation in C code{CLEAR}")
        return False
    py_solutions = output_potentials(system, ring, n1, w, True)
    if len(c_solutions) != len(py_solutions):
        print(f"{FAIL}Solution amounts differ{CLEAR}")
        print(c_solutions.keys())
        print(py_solutions.keys())
        return False
    for c_s, p_s in zip(c_solutions, py_solutions):
        if (c_s != p_s):
            print(f"{FAIL}Solutions differ!{CLEAR}", c_s, "!=", p_s)
            print(c_solutions.keys, py_solutions.keys())
            return False
    return True

def test_c_solve(sys_tuple):
    print(f"\n{ITER}== {i} =={CLEAR}")
    system, n, m, ring, known_sol = sys_tuple
    print("Verify known solution:", [f(*convert(known_sol, n)) for f in system]) 
    print("Known solution:", known_sol)
    print("System:", system, ring)
    sl_system = bitslice(system, ring.gens())
    # print(sl_system)
    c_sol = c_solve(sl_system, n, m)
    py_sol = solve(system, ring)
    print("Solution found:", c_sol)
    if c_sol == None:
        print(f"{FAIL}Found No solution{CLEAR}")
        return False
    elif (py_sol == c_sol):
        print(f"{SUCCESS}Found identical solution to sage code.{CLEAR}")
        return True
    else:
        out = [f(*c_sol) for f in system]
        print("Verify found solution:", out)
        if not all(val == 0 for val in out):
            print(f"{FAIL}Solution not valid{CLEAR}")
            return False
        print(f"{WARNING}Solution verified, but not identical to sage solution{CLEAR}")
        print(f"{py_sol} != {c_sol}")
        return False

def test_u_values(sys_tuple):
    system, n, m, ring, _ = sys_tuple
    n1 = int(ceil(n/(5.4))) 
    # print(type(n1), type(n))
    F_tilde = product((GF(2)(1) + f) for f in system)
    w = F_tilde.degree() - n1
    V, ZV = compute_u_values(system, ring, n1, w)
    for y in range(2^(n - n1)):
        print("Computing s0 sum...")
        s0 = sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1)) for z_hat in range(2^n1))
        print("Finished sum")
        if V[y] != s0:
            print(f"Error found in V[{y}]\n\t{V[y]}\n\t{s0}")
            return False
        for i in range(n1):
            print(f"Computing s{i} sum...")
            si = sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1 - 1)[:i], 0, *convert(z_hat, n1 - 1)[i:]) for z_hat in range(2^(n1 - 1)))
            print("Finished sum")
            if (ZV[i][y]) != si:
                print(f"Error found in ZV[{i}][{y}]\n\t{(ZV[i][y])}\n\t{si}")
                return False
    return True

def test_dinur_output_sol(sys_tuple):
    system, n, m, ring, _ = sys_tuple
    n1 = int(ceil(n/(5.4))) 
    F_tilde = product((GF(2)(1) + f) for f in system)
    w = F_tilde.degree() - n1
    potentials = output_potentials(system, ring, n1, w, True) 
    for y in range(2^(n - n1)):
        s0 = sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1)) for z_hat in range(2^n1))
        if s0 == GF(2)(1):
            if y not in potentials:
                print(f"Error in out: {y} is not a key")
                return False
            if potentials[y][0] != GF(2)(1):
                print(f"Error in out[{y}][0]: Not equal to 1 when sum evaluates to 1")
                return False
            for i in range(1, n1 + 1):
                si = sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1 - 1)[:i - 1], 0, *convert(z_hat, n1 - 1)[i - 1:]) for z_hat in range(2^(n1 - 1)))
                if potentials[y][i] == si:
                    print(f"Error in out[{y}][{i}]: Value NOT inverse of the evaluation of U_{i}({y}) = {si} = {potentials[y][i]}")
                    return False
    return True

def test_run_solve(sys_tuple):
    print(f"\n{ITER}== {i} =={CLEAR}")
    system, n, _, ring, known_sol = sys_tuple
    print("Verify known solution:", [f(*convert(known_sol, n)) for f in system]) 
    print("Known solution:", known_sol)
    print("System:", system, ring)
    sol = solve(system, ring)
    print("Solution found:", sol)
    if sol == None:
        print(f"{FAIL}Found No solution{CLEAR}")
        return False
    else:
        out = [f(*sol) for f in system]
        print("Verify found solution:", out)
        if not all(val == 0 for val in out):
            print(f"{FAIL}Solution not valid{CLEAR}")
            return False
        else:
            print(f"{SUCCESS}Solution verified{CLEAR}")
            return True

def compute_u_values(system, ring, n1, w):
    global _time_bruteforce
    global _time_u_values

    n = len(ring.gens())

    _time_bruteforce -= time.time() 

    sols = bruteforce(system, ring.gens(), n1, w + 1)

    _time_bruteforce += time.time()

    V, ZV = None, None

    V = defaultdict(lambda: GF(2)(0))
    ZV = [defaultdict(lambda: GF(2)(0)) for _ in range(n1)]

    _time_u_values -= time.time()

    for s in sols:
        y, z = s[:n - n1], s[n - n1:]

        if sum(y) <= w:
            idx = index_of(y)
            V[idx] += 1
        
        for i in range(1, n1 + 1):
            if z[i - 1] == 0:
                idx = index_of(y)
                ZV[i - 1][idx] += 1

    _time_u_values += time.time()

    return V, ZV

def output_potentials(system, ring, n1, w, fes_recovery):
    global _time_mobius
    global _time_full_eval
    global _time_fetch_sol
    global _time_fes_recovery

    n = len(ring.gens())

    out = defaultdict(lambda: GF(2)(0))

    if fes_recovery:
        
        _time_fes_recovery -= time.time()
        
        evals = fes_recover(system, n, n1, w + 1, ring)
        
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
        ring_sub = GF(2)[", ".join([str(var) for var in ring.gens()[:n - n1]])]
        U = []

        V, ZV = compute_u_values(system, ring, n1, w + 1)

        _time_mobius -= time.time() 

        U.append(mob_transform(V, ring_sub.gens(), w))
        for i in range(1, n1 + 1):
            U.append(mob_transform(ZV[i - 1], ring_sub.gens(), w+1))

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
    
            tmp = mob_transform(tmp, ring_sub.gens())
    
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

    return out

def c_output_potentials(system, n, n1, w, test=False):
    c_system = (C_POLY_T * len(system))(*system)
    out = (C_VARS_T * int(2^(n - n1)))()

    out_size = ct.c_size_t(0)

    args = [Type.P(C_POLY_T), Type.U, Type.U, Type.U, Type.P(C_VARS_T), Type.P(Type.SZ)]
    res = Type.U8
    output_fun = fetch_c_func("output_potentials", args, res, test)
    error = output_fun(c_system, n, n1, w, out, ct.pointer(out_size))
    
    if error == 1:
        return None
    
    py_dict = defaultdict(lambda: GF(2)(0))
    for y_hat in range(2^(n - n1)):
        if (out[y_hat] & 1) == 1:
            py_dict[y_hat] = out[y_hat]
    
    return py_dict

def eval_system(system, sol):
    return not any(f(*sol) for f in system)

def preprocess(system, ring):
    new_sys = []
    for f in system:
        for v in ring.gens():
            if v^2 in f:
                f = f + v^2 + v
        new_sys.append(f)
    return new_sys

def gen_matrix_rank_l(l, m):
    mat = matrix([convert(rand() & ((1 << m) - 1), m) for _ in range(l)])
    while mat.rank() != l:
        print(mat)
        mat = matrix([convert(rand() & ((1 << m) - 1), m) for _ in range(l)])
    print(mat)
    return mat    

def solve(system, ring, fes_recovery=True):
    srand(int(RSEED)) # Seeding randomness in accordance with C code.

    global _time_solve_trials   
    global _time_output_potentials
    system = preprocess(system, ring)

    n = len(ring.gens())
    n1 = int(ceil(n/(5.4))) # Quadratic systems are assumed here, see page 19 of full dinur paper for explanation
    l = n1 + 1
    m = len(system)
    potentials_solutions = []
    k = 0

    while True:
        print("Commencing round", k)
        A = gen_matrix_rank_l(l, m)

        E_k = [sum(GF(2)(A[i][j]) * system[j] for j in range(m)) for i in range(l)]
        w = sum(f.degree() for f in E_k) - n1 
        _time_output_potentials -= time.time()

        curr_potential_sol = output_potentials(E_k, ring, n1, w, fes_recovery) 
        # for k in curr_potential_sol:
        #     print(list(convert(k, n - n1)) + list(curr_potential_sol[k]), end=" ")
        # print()
        _time_output_potentials += time.time()

        potentials_solutions.append(curr_potential_sol)

        _time_solve_trials -= time.time()
        for y_hat, potential_sol in curr_potential_sol.items(): # Iterate through solutions instead of all possible inputs
            for k1 in range(k):
                if all(potential_sol == potentials_solutions[k1][y_hat]):
                    sol = convert(y_hat, n - n1) + list(potential_sol[1:])
                    if eval_system(system, sol):
                        _time_solve_trials += time.time()
                        # print("Rounds:", k + 1)
                        return sol
                    break
        _time_solve_trials += time.time()
        k += 1
    return None

def c_solve(system, n, m, test=False):
    c_system = (C_POLY_T * len(system))(*system)
    sol = C_VARS_T(0)

    args = [Type.P(C_POLY_T), Type.U, Type.U, Type.P(C_VARS_T)]
    res = Type.U8
    solve_fun = fetch_c_func("solve", args, res, test)
    error = solve_fun(c_system, n, m, ct.pointer(sol))

    if (error == 1):
        return None
    
    return convert(sol.value, n)

import sys

def test_c_compute_e_k(sys_tuple):
    srand(RSEED)
    system, n, m, ring, _ = sys_tuple
    l = int(ceil(n/(5.4))) + 1
    sl_system = bitslice(system, ring.gens())
    mat = gen_matrix_rank_l(l, m)
    mat_dec = [index_of(e) for e in mat]
    # TODO: Use subprocesses to run and communicate with C code.
    print(l)
    print(n)
    print(len(sl_system))
    for e in mat_dec:
        print(e, end=" ")
    for e in sl_system:
        print(e, end=" ")
    new_sys = bitslice([sum(mat[i][j] * system[j] for j in range(m)) for i in range(l)], ring.gens())
    print(sum(f.degree() for f in [sum(mat[i][j] * system[j] for j in range(m)) for i in range(l)]))
    for e in new_sys:
        print(e, end=" ")
    return True

def test_c_eval(sys_tuple):
    # TODO: Use subprocesses to run and communicate with C code.
    
    system, n, m, ring, _ = sys_tuple
   
    sl_sys = bitslice(system, ring.gens())

    print(len(sl_sys))
    print(n)
    for i in sl_sys:
        print(i, end=" ")

    for i in [index_of([int(f(*convert(i, n))) for f in system]) for i in range(2^n)]:
        print(i, end=" ")
    return True

def test_c_gen_matrix(sys_tuple):
    # TODO: Use subprocesses to run and communicate with C code.
    _, n, m, _, _ = sys_tuple 
    l = int(ceil(n/(5.4))) + 1
    print(l)
    print(m)
    mat = gen_matrix_rank_l(l, m)
    for i in mat:
        print(index_of(i), end=" ")
    return True

def main():
    rounds = 1
    # test_run_solve(sys_tuple)
    # print("Bruteforce time:", _time_bruteforce/rounds)
    # print("U value computation time:", _time_u_values/rounds)
    # print("Time for FES interpolation:", _time_fes_recovery/rounds)
    # print("Möbius transform time", _time_mobius/rounds)
    # print("Full evaluation time", _time_full_eval/rounds)
    # print("Fetch solutions time", _time_fetch_sol/rounds)
    # print("Output potentials time", _time_output_potentials/rounds)
    # print("Solve trials time:", _time_solve_trials/rounds)
    # print("Solve time", _time_solve/rounds)

if __name__ == "__main__":
    main()