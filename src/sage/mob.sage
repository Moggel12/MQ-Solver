from random import randint
from itertools import product as prod
from utils import index_of
from math import comb as c

from collections import defaultdict

def next_hw(v):
    t = (v | (v - 1)) + 1
    return t | ((((t & -t) // (v & -v)) >> 1) - 1)

# def get_sub_patterns(bits, table, indices):
#     if indices in table:
#         return table[indices]
#     new_bits = [bits ^ (1 << indices[0]), bits]
#     if len(indices) == 1:
#         table[indices] = new_bits
#         return new_bits 
#     new_indices = indices[1:]
#     table[new_indices] = get_sub_patterns(new_bits[0], table, new_indices) + get_sub_patterns(new_bits[1], table, new_indices)
#     return table[new_indices]

def get_sub_patterns(bits, indices):
    new_bits = [bits ^ (1 << indices[0]), bits]
    if len(indices) == 1:
        return new_bits
    new_indices = indices[1:]
    return get_sub_patterns(new_bits[0], new_indices) + get_sub_patterns(new_bits[1], new_indices)

def sparse_mob_transform(sol, vars, max_degree):
    n = len(vars)
    new_sol = defaultdict(lambda: GF(2)(0))
    for d in range(1, max_degree + 1):
        inp = (2 << d) - 1 # Initial bits with HW == d
        for _ in range(c(n, d) - 1): # Go through all bitstrings of HW == d
            indices = tuple(i for i in range(n) if (1 & (inp >> i)) == 1) # Compute indices of current bitstring
            for inp_reduced in get_sub_patterns(inp, indices): # Go through all sub_patterns to compute t_I
                 new_sol[inp] += sol[inp_reduced]
            inp = next_hw(inp) # Get next bitstring of HW == d
    return mob_transform(new_sol, vars)

# def _descent(A, B, C, n, d, k, b):
#     if k == 0:
#         pass
#     else:
#         pass
#         _descent(B, A, C, n - 1, k - 1, b)

def space_efficient_mob_transform(A, n, d, k):
    B = np.full(len(A), GF(2)(0))
    C = np.full(2^(n-k), GF(2)(0))
    for b in prod([0,1], repeat=k):
        _descent(A, B, C, n, d, k, b)

def mob_transform(sol, vars):
    return _f_expand(0, vars, sol)

def _f_expand(lvl, vars, sol):
    alt = lambda val : [val if lvl == idx else e for idx, e in enumerate(vars)]
    if lvl == len(vars):
        idx = index_of(vars)
        return sol[idx]
    else:
        f1 = (_f_expand(lvl + 1, alt(0), sol) + _f_expand(lvl + 1, alt(1), sol))
        f2 = _f_expand(lvl + 1, alt(0), sol)
        return vars[lvl] * f1 + f2

def mob_inv(poly):
    sol = dict()
    for idx, bits in enumerate(product([GF(2)(0), GF(2)(1)], repeat=len(poly.args()))):
        sol[idx] = poly(*(bits[::-1]))
    return sol

def random_poly_GF2(R, degree):
    return R(GF(2)[R.gens()].random_element(degree=degree))

def test_mob():
    for _ in range(10):
        num_vars = randint(2,10)
        R = GF(2)[", ".join(["x" + str(i) for i in range(num_vars)])]
        p = random_poly_GF2(R, 2)
        sol = mob_inv(p)
        p_anf = mob_transform(sol, R.gens())
        p_anf_sol = mob_inv(p_anf)
        if sol != p_anf_sol:
            print(f"p{p.args()} =", p)
            print(f"p_anf{p_anf.args()} =", p_anf)
            print(sol)
            print(p_anf_sol)
    print("Testing done")

def test_mob2():
    sol = {i: 1 for i in range(8)}
    n = 3
    R.<x0, x1, x2> = GF(2)[]
    p = mob_transform(sol, R.gens(), 2)
    print(p)

def main():
    # test_mob()
    test_mob2()

if __name__ == "__main__":
    main()
