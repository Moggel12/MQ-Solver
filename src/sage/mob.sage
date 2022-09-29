from random import randint
from itertools import product

def mob_transform(sol, vars):
    return _f_expand(0, vars, sol)

def _f_expand(lvl, vars, sol):
    alt = lambda val : [val if lvl == idx else e for idx, e in enumerate(vars)]
    if lvl == len(vars):
        idx = sum(bit << pos for (pos, bit) in enumerate(reversed(vars)))
        return sol[idx]
    else:
        f1 = (_f_expand(lvl + 1, alt(0), sol) + _f_expand(lvl + 1, alt(1), sol))
        f2 = _f_expand(lvl + 1, alt(0), sol)
        return vars[lvl] * f1 + f2

def mob_inv(poly):
    sol = dict()
    for idx, bits in enumerate(product([GF(2)(0), GF(2)(1)], repeat=len(poly.args()))):
        sol[idx] = poly(*bits)
    return sol

def random_poly_GF2(R, degree):
    return R(GF(2)[R.gens()].random_element(degree=degree))

def test_mob():
    for _ in range(100):
        num_vars = randint(2,10)
        R = GF(2)[", ".join(["x" + str(i) for i in range(num_vars)])]
        p = random_poly_GF2(R, 2)
#    print(f"p{p.args()} =", p)
        sol = mob_inv(p)
#    print(sol)
        p_anf = mob_transform(sol, R.gens())
#    print(f"p_anf{p_anf.args()} =", p_anf)
        p_anf_sol = mob_inv(p_anf)
#    print("Solutions equal:", sol == p_anf_sol)
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
