from dataclasses import dataclass

from utils import convert, bitslice, fetch_c_func
from itertools import combinations

@dataclass
class State:
    i: int
    y: int
    d1: list
    d2: list
    prefix: list

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
    n = len(vars)
    sliced = bitslice(system, vars)
    s = None
    for i in range(2^(n - n1)):
        if hamming_weight(i) > d:
            continue
        prefix = [pos for pos, b in enumerate(reversed(bin(i)[2:])) if b == "1"]
        s = update(s, sliced, n, n1, prefix)
        sub_sol = fes_eval(sliced, n, n1, prefix, s)
        solutions += [convert(sol, n) for sol in sub_sol]
    return solutions



def c_fes_eval_parities():
    c_func = fetch_c_func("test")
    c_func()
    print("####")
    system = [ 5, 5, 25, 1, 7, 12, 17, 1, 28, 29, 2, 16, 21, 15, 21, 30]
    print(system)
    prefix = [0]
    n = 5
    n1 = 3
    s = None
    s = update(s, system, n, n1, prefix)
    for s in fes_eval(system, n, n1, prefix, s, False):
        print("Solution", s)
    
if __name__ == "__main__":
    c_fes_eval_parities()
