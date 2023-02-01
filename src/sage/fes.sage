from dataclasses import dataclass

from utils import *
from itertools import combinations
import ctypes as ct

@dataclass
class State:
    i: int
    y: int
    d1: list
    d2: list
    prefix: list

class C_State(ct.Structure):
    _fields_ = [("i", ct.c_uint8),
                ("y", ct.c_uint8),
                ("d1", ct.POINTER(ct.c_uint8)),
                ("d2", ct.POINTER(ct.c_uint8)),
                ("prefix", ct.POINTER(ct.c_uint8))]

def hamming_weight(x):
    x -= (x >> 1) & 0x5555555555555555
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333)
    x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f
    return ((x * 0x0101010101010101) & 0xffffffffffffffff ) >> 56

def preprocess(f):
    for v in f.parent.gens():
        if v^2 in f.monomials():
            f = f + v^2 + v
    return f

def lex_idx(i, j, n):
    i, j = (i, j) if i < j else (j, i)
    return n + sum((n - k) for k in range(1, i + 2)) - (n - j - 1)

# Get index position of first bit set (if any).
def bit1(x):
    x = x&-x

    x = int(x).bit_length()

    return None if x == 0 else x-1

# Get index position of second bit set (if any).
def bit2(x):
    # Remove first set bit and get position of the remaining 'first' bit.
    return bit1(x ^^ (x&-x))

def init(f, n, n1, prefix):
    s = State(i = 0, y = f[0], d1=[0]*n1, d2=[[0]*n1 for _ in range(n1)], prefix=prefix)

    for k in range(n1):
        for j in range(k):
            s.d2[k][j] = f[lex_idx(j + (n - n1), k + (n - n1), n)] # Alter this function in old FES

    s.d1[0] = f[1 + n - n1]
    for k in range(1, n1):
        s.d1[k] = s.d2[k][k-1] ^^ f[1 + k + (n - n1)]

    # add pre-evaluation for prefix variables
    for idx in prefix:
        for k in range(n1):
            s.d1[k] ^^= f[lex_idx(idx, k + (n - n1), n)] # Alter this function in old FES
            
        s.y ^^= f[idx + 1]

    for i, j in combinations(prefix, 2):
        idx = lex_idx(i, j, n) # Alter this function in old FES
        s.y ^^= f[idx] # Index into this here

    return s

def update(s, f, n, n1, prefix):
    if s == None:
        return init(f, n, n1, prefix)

    off = [v for v in s.prefix if v not in prefix]
    on = [v for v in prefix if v not in s.prefix]

    # turn old variables off
    for idx in off:
        for k in range(n1):
            s.d1[k] ^^= f[lex_idx(idx, k + (n - n1), n)] 

        s.y ^^= f[idx + 1]

    for i in off:
        for j in [v for v in s.prefix if v not in off]:
            s.y ^^= f[lex_idx(i, j, n)]

    for i, j in combinations(off, 2):
        s.y ^^= f[lex_idx(i, j, n)]

    # turn new variables on
    for idx in on:
        for k in range(n1):
            s.d1[k] ^^= f[lex_idx(idx, k + (n - n1), n)]

        s.y ^^= f[idx + 1]

    for i in on:
        for j in [v for v in prefix if v not in on]:
            s.y ^^= f[lex_idx(i, j, n)]
            # s.y += f.monomial_coefficient(X[i]*X[j]) # <---- Remove when ready

    for i, j in combinations(on, 2):
        s.y ^^= f[lex_idx(i, j, n)]
        # s.y += f.monomial_coefficient(X[i]*X[j]) # <---- Remove when ready

    s.prefix = prefix

    return s

def step(s):
    s.i = s.i + 1
    k1 = bit1(s.i)
    k2 = bit2(s.i)

    if k2:
        s.d1[k1] ^^= s.d2[k2][k1]

    s.y ^^= s.d1[k1]

def fes_eval(f, n, n1 = None, prefix=[], s = None, compute_parity=False):

    if n1 == None:
        n1 = n

    if compute_parity:
        parities = 0
    else:
        res = []

    if s == None:
        s = init(f, n, n1, prefix)
    
    pre_x = sum([1<<i for i in prefix])

    if s.y == 0:
        if compute_parity:
            parities ^^= 2^(n1 + 1) - 1
        else:
            res.append(((s.i ^^ (s.i >> 1)) << (n-n1)) | pre_x)

    while s.i < 2^n1 - 1:
        step(s)

        if s.y == 0:
            if compute_parity:
                parities ^^= 1
                z = (s.i ^^ (s.i >> 1))
                for pos in range(n1):
                    if z & (1 << pos) == 0:
                        parities ^^= (1 << (pos + 1))
            else:
                res.append(((s.i ^^ (s.i >> 1)) << (n-n1)) | pre_x)

    # reset s.i to initial state
    s.i = 0
    for i in range(n1-1):
        s.d1[i] ^^= s.d2[n1-1][i]
    s.y ^^= s.d1[n1-1] ^^ s.d2[n1-1][n1-2]

    if compute_parity:
        return parities 
    return res

def bruteforce(system, vars, n1, d):
    solutions = []
    # n = len(vars)
    n = vars
    # sliced = bitslice(system, vars)
    s = None
    for i in range(2^(n - n1)):
        if hamming_weight(i) > d:
            continue
        prefix = [pos for pos, b in enumerate(reversed(bin(i)[2:])) if b == "1"]
        s = update(s, system, n, n1, prefix)
        sub_sol = fes_eval(system, n, n1, prefix, s)
        solutions += [convert(sol, n) for sol in sub_sol]
    return solutions


def c_bruteforce(system, n, n1, d):
    print("Test")
    solutions = []
    c_system = (ct.c_uint8 * len(system))(*system)
    c_solutions = (ct.c_uint8 * int((1 << n)))(0)

    args = [ct.POINTER(ct.c_uint8), ct.c_uint, ct.c_uint, ct.c_uint, ct.POINTER(ct.c_uint8)]
    res = ct.c_uint
    brute = fetch_c_func("bruteforce", args, res)
    sol_amount = brute(c_system, n, n1, d, c_solutions)

    py_list = []
    for i in range(sol_amount):
        py_list.append(int(c_solutions[i]))

    return sol_amount, py_list 

def c_fes_eval_test(trials, m=5, n=5):

    for _ in range(trials):
        # system, _, ring, n = random_system_with_sol(m, m, n, n)
        # n1 = randint(1, n - 1)
        # d = randint(1, n - n1)
        ring.<x0,x1,x2,x3,x4> = GF(2)[]
        system = [x0*x2 + x2*x3 + x0*x4 + x2*x4 + x0 + x2 + x3, x0*x2 + x1*x2 + x0*x3 + x1*x3 + x2*x3 + x0*x4 + x1*x4 + x2*x4 + x1, x0*x1 + x0*x2 + x1*x2 + x0*x3 + x1*x3 + x2*x3 + x0*x4 + x1*x4 + x3*x4 + x0 + x1 + x2 + x3 + x4 + 1, x0*x1 + x0*x2 + x1*x3 + x2*x3 + x1*x4 + x3*x4 + x0 + x1 + x2, x0*x2 + x0*x4 + x1*x4 + x3*x4 + x2 + x3]
        n = 5
        n1 = 1
        d = 3
        print(system)
        print(n, n1, d)
        system = bitslice(system, ring.gens())
        # print("\n=== Running C code ===\n")
        amount, c_solutions = c_bruteforce(system, n, n1, d)  
        # print("\n=== Running Python code ===\n")
        py_solutions = [index_of(sol) for sol in bruteforce(system, n, n1, d)]
        c_solutions.sort()
        py_solutions.sort()
        if amount != len(py_solutions):
            print("Amount of solutions differ")
            print(c_solutions)
            print(py_solutions)
            return
        for s_c, s_p in zip(c_solutions, py_solutions):
            if s_c != s_p:
                print("Solutions differ!", s_c, "!=", s_p)
                print(c_solutions, py_solutions)
                return
    print(f"Ran {trials} trials without errors!")
    
if __name__ == "__main__":
    c_fes_eval_test(10)
