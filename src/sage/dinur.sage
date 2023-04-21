import numpy as np
from random import randint
from math import ceil
import time
import ctypes as ct
from collections import defaultdict

from src.sage.mob_new import mob_transform
from src.sage.fes import bruteforce
from src.sage.utils import index_of, convert, random_systems, fetch_systems_interactive, CLEAR, WARNING, FAIL, SUCCESS, ITER
from src.sage.c_config import Type, srand, rand, RSEED, C_POLY_T, C_VARS_T, MAX_HISTORY, C_VECTORIZED
from src.sage.fes_rec import fes_recover
from src.sage.utils import bitslice, fetch_c_func, write_fukuoka, run_bin_test

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
    if C_VECTORIZED:
        print(f"{FAIL}This test does not support vectorized C functionality.{CLEAR}")
        return False
    system, n, m, ring, _ = sys_tuple
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
            print(f"{FAIL}Keys in solutions differ!{CLEAR}", c_s, "!=", p_s)
            print(c_solutions.keys, py_solutions.keys())
            return False
        elif any(convert(c_solutions[c_s], n1 + 1) != py_solutions[p_s]):
            print(f"{FAIL}Z bits differ!{CLEAR}", py_solutions[p_s], "!=", convert(c_solutions[c_s], n1 + 1), py_solutions[p_s])
            print(c_solutions.values())
            print(py_solutions.values())
            return False
    return True

def test_c_solve(sys_tuple):
    system, n, m, ring, known_sol = sys_tuple
    print("Verify known solution:", [f(*convert(known_sol, n)) for f in system]) 
    print("Known solution:", known_sol)
    print("System:", system, ring)
    sl_system = bitslice(system, ring.gens())
    try:
        c_sol = c_solve(sl_system, n, m)
    except:
        print(f"{FAIL}Encountered exception in C code.{CLEAR}")
        return False 
    print("Solution found:", c_sol)
    if c_sol == None:
        print(f"{FAIL}Found No solution{CLEAR}")
        return False
    else:
        out = [f(*c_sol) for f in system]
        print("Verify found solution:", out)
        if not all(val == 0 for val in out):
            print(f"{FAIL}Solution not valid{CLEAR}")
            return False
        print(f"{SUCCESS}Solution verified.{CLEAR}")
        return True

def test_sage_u_values(sys_tuple):
    system, n, m, ring, _ = sys_tuple
    n1 = int(ceil(n/(5.4))) 
    F_tilde = product((GF(2)(1) + f) for f in system)
    w = F_tilde.degree() - n1
    V, ZV = compute_u_values(system, ring, n1, w)
    for y in range(2^(n - n1)):
        s0 = sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1)) for z_hat in range(2^n1))
        if V[y] != s0:
            print(f"Error found in V[{y}]\n\t{V[y]}\n\t{s0}")
            return False
        for i in range(n1):
            si = sum(F_tilde(*convert(y, n - n1), *convert(z_hat, n1 - 1)[:i], 0, *convert(z_hat, n1 - 1)[i:]) for z_hat in range(2^(n1 - 1)))
            if (ZV[i][y]) != si:
                print(f"Error found in ZV[{i}][{y}]\n\t{(ZV[i][y])}\n\t{si}")
                return False
    return True

def test_sage_output_solutions(sys_tuple):
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

def test_sage_solve(sys_tuple):
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
                    out[y_hat][i] = GF(2)((evals[y_hat] >> i) & 1) + 1 # Bitsliced indexing

        _time_fetch_sol += time.time()
    else:
        ring_sub = GF(2)[", ".join([str(var) for var in ring.gens()[:n - n1]])]
        U = []

        V, ZV = compute_u_values(system, ring, n1, w)

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

    for idx in range(out_size.value):
        y_hat = out[idx] & (2^(n - n1) - 1)
        z_bits = (out[idx] >> (n - n1 - 1)) | 1
        py_dict[y_hat] = z_bits
    
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

def gen_matrix_rank_l(l, m, vec=False):
    mat = matrix([convert(rand() & ((1 << m) - 1), m) for _ in range(l)])
    while mat.rank() != l:
        mat = matrix([convert(rand() & ((1 << m) - 1), m) for _ in range(l)])
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
    potential_solutions = []
    k = 0

    while k < MAX_HISTORY:
        print("Commencing round", k)
        A = gen_matrix_rank_l(l, m)
        
        E_k = [sum(GF(2)(A[i][j]) * system[j] for j in range(m)) for i in range(l)]
        
        w = sum(f.degree() for f in E_k) - n1 

        _time_output_potentials -= time.time()

        curr_potential_sol = output_potentials(E_k, ring, n1, w, fes_recovery) 

        _time_output_potentials += time.time()

        potential_solutions.append(curr_potential_sol)

        _time_solve_trials -= time.time()
        for y_hat, potential_sol in curr_potential_sol.items(): # Iterate through solutions instead of all possible inputs
            if potential_sol[0] == 1:
                for k1 in range(k):
                    if all(potential_sol == potential_solutions[k1][y_hat]):
                        sol = convert(y_hat, n - n1) + list(potential_sol[1:])
                        if eval_system(system, sol):
                            #print(y_hat, index_of(list(potential_sol[1:])), index_of(sol))
                            _time_solve_trials += time.time()
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

def c_benchmark(system_tuples):
    n = system_tuples[0][1]
    m = system_tuples[0][2]
    ring = system_tuples[0][3]
    size = len(system_tuples)

    c_systems_list = (Type.P(C_POLY_T) * size)()
    for i in range(size):
        sl_sys = bitslice(system_tuples[i][0], ring.gens())
        c_systems_list[i] = (C_POLY_T * int(math.comb(n, 2) + n + 1))(*sl_sys)

    args = [Type.SZ, Type.P(Type.P(C_POLY_T)), Type.SZ, Type.SZ]
    bench_fun = fetch_c_func("e2e_benchmark", args)
    bench_fun(size, c_systems_list, Type.SZ(n), Type.SZ(m))

def test_c_compute_e_k_SAN(sys_tuple):
    srand(int(RSEED))  
    if C_VECTORIZED:
        amnt = 4
    else:
        amnt = 1
    system, n, m, ring, _ = sys_tuple
    l = int(ceil(n/(5.4))) + 1
    sl_system = bitslice(system, ring.gens())
    mat = []
    mat_dec = []
    for v in range(amnt):
        mat.append(gen_matrix_rank_l(l, m))
        mat_dec.append([index_of(e) for e in mat[v]])
    input_data = f"{0}\n"
    input_data += f"{l}\n{n}\n{len(sl_system)}\n"
    for i in range(amnt):
        for e in mat_dec[i]:
            input_data += f"{e} "
    for e in sl_system:
        input_data += f"{e} "
    new_sys = []
    for v in range(amnt):
        new_sys.append(bitslice([sum(mat[v][i][j] * system[j] for j in range(m)) for i in range(l)], ring.gens()))
        input_data += f"{sum(f.degree() for f in [sum(mat[v][i][j] * system[j] for j in range(m)) for i in range(l)])}\n"
    for v in range(amnt):
        for e in new_sys[v]:
            input_data += f"{e} "
    p = run_bin_test(input_data)
    if (not p) or (p.returncode != 0):
        return False
    return True

def test_c_eval_SAN(sys_tuple):    
    system, n, m, ring, _ = sys_tuple
    sl_sys = bitslice(system, ring.gens())
    input_data = f"{1}\n"
    input_data += f"{len(sl_sys)}\n{n}\n"
    for i in sl_sys:
        input_data += f"{i} "
    for i in [index_of([int(f(*convert(i, n))) for f in system]) for i in range(2^n)]:
        input_data += f"{i} "
    p = run_bin_test(input_data)
    if (not p) or (p.returncode != 0):
        return False
    return True

def test_c_gen_matrix_SAN(sys_tuple):
    srand(int(RSEED))
    if C_VECTORIZED:
        amnt = 4
    else:
        amnt = 1
    _, n, m, _, _ = sys_tuple 
    l = int(ceil(n/(5.4))) + 1
    input_data = f"{2}\n"
    input_data += f"{l}\n{m}\n"
    mat = []
    for _ in range(amnt):
        mat.append(gen_matrix_rank_l(l, m))
    for v in range(amnt):
        for i in mat[v]:
            input_data += f"{index_of(i)} "
    p = run_bin_test(input_data)
    if (not p) or (p.returncode != 0):
        return False
    return True

def test_c_solve_SAN(sys_tuple):
    system, n, m, ring, _ = sys_tuple
    solution = solve(system, ring)
    if solution == None:
        print(f"{FAIL}Sage code could not find a solution.{CLEAR}")
        return False
    sl_system = bitslice(system, ring.gens())
    input_data = f"{3}\n"
    input_data += f"{n}\n{m}\n{len(sl_system)}\n"
    if not C_VECTORIZED:
        input_data += f"{index_of(solution)} "
    for e in sl_system:
        input_data += f"{e} "
    p = run_bin_test(input_data)
    if (not p) or (p.returncode != 0):
        return False
    return True

def main():
    rounds = 1
    # print(test_c_output_solutions(sys_tuple))
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