import numpy as np
from itertools import combinations
import math
import re

_REGEX_FIELD = r"Galois Field : GF\(2\)"
_REGEX_NUM_VARS = r"Number of variables \(n\) : (\d+)"
_REGEX_NUM_POLYS = r"Number of equations \(m\) : (\d+)"
_REGEX_SYS = r"(([0-1] )+);"
_MQ_URL = "www.mqchallenge.org/format.html"
_WARNING = "\033[93m"
_FAIL = "\033[91m"
_CLEAR = "\033[0m"

def bitslice(f_sys, vars):
    f_sys_sliced = np.zeros(math.comb(len(vars), 2) + len(vars) + 1, dtype=int)
    for j, poly in enumerate(f_sys):
        if poly in [GF(2)(0), GF(2)(1)]:
            f_sys_sliced[0] ^^= int(poly) << j
        else:
            f_sys_sliced[0] ^^= int(poly.constant_coefficient()) << j
            i = 1
            for v in vars:
                f_sys_sliced[i] ^^= int(poly.coefficient({v_: 1 if v == v_ else 0 for v_ in vars})) << j
                i += 1
            for v1, v2 in combinations(range(len(vars)), 2):
                f_sys_sliced[i] ^^= int(poly.coefficient({vars[v1]: 1, vars[v2]: 1, **{v: 0 for v in vars if v not in [vars[v1], vars[v2]]}})) << j
                i += 1
    return f_sys_sliced

def convert(v, n):
    v = bin(v)[2:]
    return list(map(lambda i : (int(v[-i]) if i <= len(v) else 0), range(1, n + 1)))

def index_of(y_list):
    return sum(b << i for i, b in enumerate(y_list))

# Reads monomials such that all monomials of the same are read right->left (low to high using the graded reverse lex order).
# Could also yield the values needed for easy bitslicing early on, however, this will do for now.
def read_poly(poly_string, num_vars):
    poly_string = poly_string.replace(" ", "")

    ring = GF(2)[", ".join(["x" + str(i) for i in range(num_vars)])]
    X = ring.gens()
    f = GF(2)(poly_string[-1]) 

    num_quads = math.comb(num_vars + 2 - 1, 2)
    quad_i_j_idx = num_quads - 1 # Ensures quadratics are read left->right

    for i in range(num_vars - 1, -1, -1):
        linear_i_idx = i + num_quads # Indexes left->right as i=(num_vars-1)..0
        for j in range(i, -1, -1):
            if i == j:
                f += GF(2)(poly_string[quad_i_j_idx]) * X[i]
            else:
                f += GF(2)(poly_string[quad_i_j_idx]) * X[i] * X[j]
            quad_i_j_idx -= 1
        f += GF(2)(poly_string[linear_i_idx]) * X[i]
    return f

def read_contents(contents):
    system = [0, 0, []]
    if re.search(_REGEX_FIELD, contents) == None:
        print(f"{_FAIL}Polynomial system must be declared as GF(2).{_CLEAR} ({_MQ_URL})")
        return None
    num_vars = re.search(_REGEX_NUM_VARS, contents)
    if num_vars == None:
        print(f"{_FAIL}Number of variables is not declared correctly.{_CLEAR} ({_MQ_URL})")
        return None
    system[0] = int(num_vars[1])
    num_polys = re.search(_REGEX_NUM_POLYS, contents)
    if num_polys == None:
        print(f"{_FAIL}Number of polynomials is not declared correctly.{_CLEAR} ({_MQ_URL})")
        return None
    system[1] = int(num_polys[1])
    if re.search(_REGEX_SYS, contents) == None:
        print(f"{_FAIL}Polynomial declarations are not formatted correctly.{_CLEAR} ({_MQ_URL})")
        return None
    l = 0
    for p in re.finditer(_REGEX_SYS, contents):
        system[2].append(read_poly(p[1], system[0]))
        l += 1
    if l != system[1]:
        print(f"{_WARNING}Warning: File declared {system[1]} polynomials but defined {len(system[2])}{_CLEAR}")
        system[1] = l
    return system

def parse_fukuoka(file_path, print_system=False):
    print("This functionality assumes graded reverse lex order as the monomial order for all polynomials in the system.")
    with open(file_path, "r") as file:
        contents = file.read()
        system = read_contents(contents)
        if system == None:
            print("Please run parser again after the input file has been fixed")
            return
    if print_system:
        print(f"System consists of {system[1]} polynomials in {system[0]} variables.")
        print("Polynomials:")
        for poly in system[2]:
            print(poly)
    return system
