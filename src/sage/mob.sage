from random import randint
from itertools import product
from utils import index_of

def mob_transform(sol, vars):
    return _f_expand(0, vars, sol)

def _f_expand(lvl, vars, sol):
    alt = lambda val : [val if lvl == idx else e for idx, e in enumerate(vars)]
    if lvl == len(vars):
        idx = index_of(vars)
        return sol[idx]
    else:
        tmp = _f_expand(lvl + 1, alt(0), sol)
        f1 = tmp + _f_expand(lvl + 1, alt(1), sol)
        f2 = tmp
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

def main():
    test_mob()

if __name__ == "__main__":
    main()
