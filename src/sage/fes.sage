import math
<<<<<<< Updated upstream
from random import randint

def run_fes(f_sys, n):
  for f in f_sys: # Should be changed to the "clamped" version
    s = init(f, n)
    while s["i"] < 2^n:
      next(s)
      if s["y"] == 0:
        return s["y"]

def next(s):
  s["i"] += 1
  k1 = BIT1(s["i"])
  k2 = BIT2(s["i"])
  if math.log2(s["i"]).is_integer():
    s["deriv1"][k1] = (s["deriv1"][k1] + s["deriv2"][k1][k2])
  s["y"] = (s["y"] + s["deriv1"][k1])

def init(f, n):
  s = dict()
  s["i"] = 0
  s["x"] = 0
  s["y"] = f.constant_coefficient()
  s["deriv2"] = dict()
  for k in range(1, n):
    s["deriv2"][k] = dict()
    for j in range(k):
      s["deriv2"][k][j] = f.derivative(f.variables()[k], f.variables()[j])
  s["deriv1"] = dict()
  s["deriv1"][0] = f.coefficient({v: (1 if i == 0 else 0) for i, v in enumerate(f.variables())})
  for k in range(1, n):
    s["deriv1"][k] = s["deriv2"][k][k-1] + f.coefficient({v: (1 if i == k else 0) for i, v in enumerate(f.variables())})
  return s


def BIT1(i):
  for idx, c in enumerate(reversed(bin(i ^ (i >> 1))[2:])):
    if c == 1: return idx

def BIT2(i):
  return len(bin(i ^ (i >> 1))[2:]) - 1

def main():
  for _ in range(10):
    m = randint(5, 10)
    n = randint(2, 5)
    R = GF(2)[", ".join(["x" + str(i) for i in range(n)])]
    f_sys = []
    for i in range(m):
      f_sys.append(R(GF(2)[R.gens()].random_element(degree=2)))
    run_fes(f_sys, n)

if __name__ == "__main__":
  main() 
=======
# from random import randint
import numpy as np
from itertools import combinations_with_replacement as cwr
from itertools import product as prod

def run_fes(f_sys, vars):
    solutions = []
    n = len(vars)
    s = init(f_sys, vars)
    if s["y"] == 0:
        solutions.append(0)
    while s["i"] < 2^n - 1:
        s = next(s)
        if (s["y"] == 0).all():
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

def slice_int(f_sys, vars, m):
    f_sys_sliced = np.zeros(math.comb(len(vars) + 2 - 1, 2) + len(vars) + 1, dtype=int)
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
            v1 = values[i]
            f_sys_eval[0] = f_sys_eval[0] ^^ (v0 * v1 * f_sys[lex_idx(i, j, n) + n + 1]) # Add evaluated linear terms
        for j in range(N, n):
            f_sys_eval[j - N + 1] = f_sys_eval[j - N + 1] ^^ f_sys[lex_idx(i, j, n) + n + 1]
    f_sys_eval = np.append(f_sys_eval, f_sys[lex_idx(N, N, n) + n + 1:]) # Append square terms
    return f_sys_eval
            
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

def convert(v, n):
    v = bin(v)[2:]
    return tuple(map(lambda i : (int(v[-i]) if i <= len(v) else 0), range(1, n + 1)))

def test_solutions(f_sys, sol, R):
    for s in range(2^len(R.gens())):
        faulted_sol = False
        for f in f_sys:
            args = convert(s, len(R.gens()))
            res = f(*args) 
            faulted_sol = faulted_sol or res
            if res and (s in sol):
                print("=== Case 1 ===")
                print(R.gens())
                print(f)
                print(args)
                print(res)
                print(s)
                print("===")
                return False
        if (not faulted_sol) and (s not in sol):
            print("=== Case 2 ===")
            print(R.gens())
            print(convert(s, len(R.gens())))
            print(s)
            print([f(*convert(s, len(R.gens()))) for f in f_sys])
            print("===")
            return False
    return True

def bruteforce(system, R, n, d):
    return [convert(s, len(R.gens())) for s in run_fes(system, R.gens())]

def preprocess(f_sys, vars):
    for i, f in enumerate(f_sys):
        for x in f.variables():
            coeff = f.coefficient({x:2})
            if coeff == 1:
                f_sys[i] += x
    return f_sys

def main():
    res = True
    for i in range(100):
        f_sys = []
        m = randint(5, 10)
        n = randint(2, 10)
        R = GF(2)[", ".join(["x" + str(i) for i in range(n)])]
        for _ in range(m):
            f = R(GF(2)[R.gens()].random_element(degree=2))
            f_sys.append(f)
        f_sys_prep = preprocess(f_sys.copy(), R.gens())
        m = len(f_sys)
        f_sys_sl = slice_int(f_sys_prep, R.gens(), m)
        sol = run_fes(f_sys_sl, R.gens())
        res = test_solutions(f_sys, sol, R)
        if not res:
            print(f_sys_prep)
            print(f_sys)
            print(f_sys_sl)
            print(sol)
            print(res)
            break
    if res:
        print("No errors found for 100 trials")

if __name__ == "__main__":
    #main()
>>>>>>> Stashed changes
