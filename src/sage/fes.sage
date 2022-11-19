import math
from random import randint
from utils import convert
import numpy as np
from itertools import combinations_with_replacement as cwr
from monotonic_gray import monotonic_bounded

def run_fes(f_sys, vars):
    solutions = []
    n = len(vars)
    # f_sys = preprocess(f_sys, n)
    s = init(f_sys, vars)
    if s["y"] == 0:
        solutions.append(0)
    while s["i"] < 2^n - 1:
        s = next(s)
        if s["y"] == 0:
            solutions.append(gray_code(s["i"]))
    return solutions 

def next(s):
    s["i"] += 1
    k1 = bit1(s["i"])
    k2 = bit2(s["i"])
    if k2 > -1:
        s["d1"][k1] = (s["d1"][k1] ^^ s["d2"][k1,k2])
    s["y"] = (s["y"] ^^ s["d1"][k1])
    return s

def init(f, vars):
    n = len(vars)
    s = dict()
    s["i"] = 0 
    s["y"] = f[0] # Updated
    s["d2"] = np.zeros((n, n), dtype=int) # Updated
    id = [(k, j) for (k,j) in cwr(range(n), 2)]
    for i, (k, j) in enumerate(id):
        if k == j: continue
        s["d2"][k,j] = f[i + len(vars) + 1] # Updated: High probability of error
    s["d1"] = np.zeros(n, dtype=int) # Updated
    s["d1"][0] = f[1] # Updated
    for k in range(1, n):
        s["d1"][k] = s["d2"][k-1,k] ^^ f[k + 1]
    return s

def bitslice(f_sys, vars):
    f_sys_sliced = np.zeros(math.comb(len(vars) + 1, 2) + len(vars) + 1, dtype=int)
    for j, poly in enumerate(f_sys):
        f_sys_sliced[0] += int(poly.constant_coefficient()) << j
        i = 1
        for v in vars:
            f_sys_sliced[i] += int(poly.coefficient({v_: 1 if v == v_ else 0 for v_ in vars})) << j
            i += 1
        for v1, v2 in cwr(range(len(vars)), 2):
            if v1 == v2:
                f_sys_sliced[i] += int(poly.coefficient({vars[v1]: 2, **{v: 0 for v in vars if v != vars[v1]}})) << j
            else:
                f_sys_sliced[i] += int(poly.coefficient({vars[v1]: 1, vars[v2]: 1, **{v: 0 for v in vars if v != vars[v1] and v != vars[v2]}})) << j
            i += 1
    return f_sys_sliced

def partial_eval(f_sys, values, n):
    N = len(values)
    f_sys_eval = [f_sys[0], *f_sys[(N + 1):(n + 1)]] # Append constants
    for i, v0 in enumerate(values):
        f_sys_eval[0] = f_sys_eval[0] ^^ (v0 * f_sys[i + 1])
        for j in range(i, N):
            v1 = values[j]
            f_sys_eval[0] = f_sys_eval[0] ^^ (v0 * v1 * f_sys[lex_idx(i, j, n) + n + 1]) # Add evaluated linear terms
        for j in range(N, n):
            f_sys_eval[j - N + 1] = f_sys_eval[j - N + 1] ^^ (v0 * f_sys[lex_idx(i, j, n) + n + 1])
    f_sys_eval = np.append(f_sys_eval, f_sys[lex_idx(N, N, n) + n + 1:]) # Append square terms
    return f_sys_eval

def bruteforce(system, R, n1, d):
    solutions = []
    # m = len(system)
    n = len(R.gens())
    sliced = bitslice(system, R.gens())
    for i in range(2^(n - n1)):
        if hamming_weight(i) > d:
            continue
        pe_sliced = partial_eval(sliced, convert(i, n - n1), n)
        sub_sol = run_fes(pe_sliced, R.gens()[(n - n1):])
        if sub_sol:
            solutions += [convert(i, n - n1) + convert(s, n1) for s in sub_sol]
    return solutions

def hamming_weight(x):
    x -= (x >> 1) & 0x5555555555555555
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
    x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f
    return ((x * 0x0101010101010101) & 0xffffffffffffffff ) >> 56

def preprocess(f_sys, n): 
    for i in range(n):
        f_sys[i + 1] = f_sys[i + 1] ^^ f_sys[lex_idx(i, i, n) + n + 1]
    return f_sys

def lex_idx(i, j, n):
    return sum((n - k) for k in range(i + 1)) - (n - j)

def bit1(i):
    return int(math.log2(i & (-i)))

def bit2(i):
    if math.log2(i).is_integer(): return -1
    i &= (i - 1)
    i &= (-1)
    return int(bit1(i))

def gray_code(i): return i ^^ (i >> 1)

def subspace_gen(n, n1, w):
    for g in monotonic_bounded(n - n1, w + 1): 
        for i in range(2^n1):
            yield g + convert(i, n1)

def test_solutions(f_sys, sol, R):
    for s in range(2^len(R.gens())):
        faulted_sol = False
        for f in f_sys:
            args = convert(s, len(R.gens()))
            res = f(*args) 
            faulted_sol = faulted_sol or res
            if res and (s in sol):
                print(
                    "=== Case 1 ===\n"
                    f"Variables: {R.gens()}\n"
                    f"Polynomial: {f}\n"
                    f"Input: {args}\n"
                    f"Actual sol: {res}\n"
                    f"Decimal input: {s}\n"
                    "==="
                )
                return False
        if (not faulted_sol) and (s not in sol):
            print(
                "=== Case 2 ===\n"
                f"Variables: {R.gens()}\n"
                f"Input: {convert(s, len(R.gens()))}\n"
                f"Decimal input: {s}\n"
                f"Actual solution: {[f(*convert(s, len(R.gens()))) for f in f_sys]}\n"
                "==="
            )
            return False
    return True

def test_fes(trials):
    res = True
    for _ in range(trials):
        f_sys = []
        m = randint(5, 10)
        n = randint(2, 10)
        R = GF(2)[", ".join(["x" + str(i) for i in range(n)])]
        for _ in range(m):
            f = R(GF(2)[R.gens()].random_element(degree=2))
            f_sys.append(f)
        m = len(f_sys)
        f_sys_sl = bitslice(f_sys, R.gens())
        sol = run_fes(f_sys_sl, R.gens())
        res = test_solutions(f_sys, sol, R)
        if not res:
            print(
                # f"{f_sys_prep}\n"
                f"{f_sys}\n"
                f"{f_sys_sl}\n"
                f"{sol}\n"
                f"{res}"
            )
            break
    if res:
        print(f"No errors found for {trials} trials")

def main():
    test_fes(100)

if __name__ == "__main__":
    main()
