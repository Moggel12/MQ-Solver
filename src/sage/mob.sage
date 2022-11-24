from random import randint
from itertools import product
from utils import index_of

def mob_transform(sol, vars, degree=Infinity):
  f = _f_expand(0, vars, sol, degree)

  # How can we avoid the high terms appearing during _f_expand?
  for m in f.monomials():
    if m.degree() > degree:
      f += m

  return f

def _f_expand(lvl, vars, sol, degree):
    alt = lambda val : [val if lvl == idx else e for idx, e in enumerate(vars)]

    if lvl != len(vars):
      if str(vars[:lvl]).count('1') > degree:
        return 0

    if lvl == len(vars):
        idx = index_of(vars)
        return sol[idx]
    else:
        tmp0 = _f_expand(lvl + 1, alt(0), sol, degree)
        tmp1 = _f_expand(lvl + 1, alt(1), sol, degree)

        f1 = tmp0 + tmp1
        f2 = tmp0

        return vars[lvl] * f1 + f2

def mob_inv(poly):
    sol = dict()
    for idx, bits in enumerate(product([GF(2)(0), GF(2)(1)], repeat=len(poly.args()))):
        sol[idx] = poly(*(bits[::-1]))
    return sol

def random_poly_GF2(R, degree):
    return R(GF(2)[R.gens()].random_element(degree=degree, terms=Infinity))

def test_mob():
    for _ in range(10):
        num_vars = randint(2,10)
        degree = 2

        R = GF(2)[", ".join(["x" + str(i) for i in range(num_vars)])]

        p = random_poly_GF2(R, degree)

        # degree-2 only
        for v in R.gens():
          if v^2 in p.monomials():
            p = p + v^2 + v

        sol = mob_inv(p)
        p_anf = mob_transform(sol, R.gens(), degree)
        p_anf_sol = mob_inv(p_anf)
        if p != p_anf:
            print(f"p{p.args()} =", p)
            print(f"p_anf{p_anf.args()} =", p_anf)
            print(sol)
            print(p_anf_sol)
            break
        else:
          print("ok!")
    print("Testing done")

def main():
    test_mob()

if __name__ == "__main__":
    main()
