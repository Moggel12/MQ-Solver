

# This file was *autogenerated* from the file mob_new.sage
from sage.all_cmdline import *   # import sage library

_sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_2 = Integer(2); _sage_const_10 = Integer(10)
from random import randint
from itertools import product as prod 
from math import comb as c
from collections import defaultdict

from src.sage.utils import index_of

def next_hw(v):
    t = (v | (v - _sage_const_1 )) + _sage_const_1 
    return t | ((((t & -t) // (v & -v)) >> _sage_const_1 ) - _sage_const_1 )

def get_sub_patterns(bits, indices):
    new_bits = [bits ^ (_sage_const_1  << indices[_sage_const_0 ]), bits]
    if len(indices) == _sage_const_1 :
        return new_bits
    new_indices = indices[_sage_const_1 :]
    return get_sub_patterns(new_bits[_sage_const_0 ], new_indices) + get_sub_patterns(new_bits[_sage_const_1 ], new_indices)

def stringify_bits(i, n):
    b_i = bin(i)[_sage_const_2 :]
    leading_zeros = n - len(b_i)
    return "0"*leading_zeros + b_i

def sparse_mob_transform(sol, vars, max_degree):
    n = len(vars)
    new_sol = defaultdict(lambda: GF(_sage_const_2 )(_sage_const_0 ))
    for d in range(_sage_const_1 , max_degree + _sage_const_1 ):
        print("d:", d)
        inp = (_sage_const_1  << d) - _sage_const_1  # Initial bits with HW == d
        for _ in range(c(n, d)): # Go through all bitstrings of HW == d
            indices = tuple(i for i in range(n) if (_sage_const_1  & (inp >> i)) == _sage_const_1 )
            for inp_reduced in get_sub_patterns(inp, indices): # Go through all sub_patterns to compute t_I
                new_sol[inp] += sol[inp_reduced]
            inp = next_hw(inp) # Get next bitstring of HW == d
    print(new_sol)
    return mob_transform(new_sol, vars)

def mob_transform(sol, vars, degree=Infinity):
  f = _f_expand(_sage_const_0 , vars, sol, _sage_const_0 , degree)

  # How can we avoid the high terms appearing during _f_expand?
  for m in f.monomials():
    if m.degree() > degree:
      f += m

  return f

def _f_expand(lvl, vars, sol, weight, degree):
    alt = lambda val : [val if lvl == idx else e for idx, e in enumerate(vars)]

    if weight > degree:
      return _sage_const_0 

    if lvl == len(vars):
        idx = index_of(vars)
        return sol[idx]
    else:
        tmp0 = _f_expand(lvl + _sage_const_1 , alt(_sage_const_0 ), sol, weight, degree)
        tmp1 = _f_expand(lvl + _sage_const_1 , alt(_sage_const_1 ), sol, weight+_sage_const_1 , degree)

        f1 = tmp0 + tmp1
        f2 = tmp0

        return vars[lvl] * f1 + f2

def mob_inv(poly):
    sol = dict()
    for idx, bits in enumerate(prod([GF(_sage_const_2 )(_sage_const_0 ), GF(_sage_const_2 )(_sage_const_1 )], repeat=len(poly.args()))):
        sol[idx] = poly(*(bits[::-_sage_const_1 ]))
    return sol

def random_poly_GF2(R, degree):
    return R(GF(_sage_const_2 )[R.gens()].random_element(degree=degree, terms=Infinity))

def test_mob():
    for _ in range(_sage_const_10 ):
        num_vars = randint(_sage_const_2 ,_sage_const_10 )
        degree = _sage_const_2 

        R = GF(_sage_const_2 )[", ".join(["x" + str(i) for i in range(num_vars)])]
        print(R)

        p = random_poly_GF2(R, degree)
        print(p)
        # degree-2 only
        for v in R.gens():
          if v**_sage_const_2  in p.monomials():
            p = p + v**_sage_const_2  + v

        sol = mob_inv(p)
        for k, v in sol.items():
            if sum(int(i) for i in bin(k)[_sage_const_2 :]) > degree:
                sol[k] = GF(_sage_const_2 )(_sage_const_0 )

        p_anf_2 = mob_transform(sol, R.gens(), degree)
        p_anf_1 = sparse_mob_transform(sol, R.gens(), degree)
        p_anf_sol_1 = mob_inv(p_anf_1)
        p_anf_sol_2 = mob_inv(p_anf_2)
        if (p != p_anf_2) or (p != p_anf_1):
            print(f"p{p.args()} =", p)
            print(f"p_anf{p_anf_1.args()} =", p_anf_1)
            print(f"p_anf{p_anf_2.args()} =", p_anf_2)
            print(sol)
            print(p_anf_sol_1)
            print(p_anf_sol_2)
            break
        else:
          print("ok!")
    print("Testing done")

def main():
    test_mob()

if __name__ == "__main__":
    main()

