import numpy as np
from itertools import combinations
import math
import re
import ctypes
import subprocess
import os

from src.sage.c_config import TEST_BIN_AVAILABLE

_REGEX_FIELD = r"Galois Field : GF\(2\)"
_REGEX_NUM_VARS = r"Number of variables \(n\) : (\d+)"
_REGEX_NUM_POLYS = r"Number of equations \(m\) : (\d+)"
_REGEX_SYS = r"(([0-1] )+);"
_REGEX_SOL = r"Known solution: (\d+)"

_MQ_URL = "www.mqchallenge.org/format.html"

WARNING = "\033[93m"
FAIL = "\033[91m"
CLEAR = "\033[0m"
SUCCESS = "\x1b[32m"
ITER = "\x1b[36m"

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
    return sum(int(b) << i for i, b in enumerate(y_list))

def fetch_c_func(func_name, args=None, res=None, test=False):
    sage_dir = os.path.dirname(__file__)
    if test: 
        libc = ctypes.CDLL(os.path.join(sage_dir, "../../bin/test.so"))
    else: 
        libc = ctypes.CDLL(os.path.join(sage_dir, "../../bin/mq.so"))
    c_func = eval(f"libc.{func_name}")
    if args: c_func.argtypes = args
    if res: c_func.restype = res
    return c_func

def run_bin_test(input_data):
    file_test = os.path.join(os.path.dirname(__file__), "../../bin/test")
    if not TEST_BIN_AVAILABLE:
        print(f"{FAIL}Test binary is not available.{CLEAR}")
        return None
    return subprocess.run([file_test], input=input_data, text=True)

def create_poly_str(poly, ring):
    X = ring.gens()
    n = len(ring.gens())

    poly_str = ""

    for i in range(n):
        for j in range(i + 1):
            poly_str += f"{poly.monomial_coefficient(X[j]*X[i])} "
    for i in range(n):
        poly_str += f"{poly.monomial_coefficient(X[i])} "
    poly_str += f"{poly.constant_coefficient()} ;\n"

    return poly_str

def write_fukuoka(filename, system, n, m, ring, sol=None):
    if sol and (type(sol) != int):
        print(f"{FAIL}Solution must be an int representing the decimal representation of the c input. Got {type(sol)}{CLEAR}")
        return
    with open(filename, "w") as f:
        f.write("Galois Field : GF(2)\n")
        f.write(f"Number of variables (n) : {n}\n")
        f.write(f"Number of equations (m) : {m}\n")
        if sol:
            f.write(f"Known solution: {sol}\n")
        f.write(f"\n*********************\n")
        for poly in system:
            f.write(create_poly_str(poly, ring))

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
        print(f"{FAIL}Polynomial system must be declared as GF(2).{CLEAR} ({_MQ_URL})")
        return None
    num_vars = re.search(_REGEX_NUM_VARS, contents)
    if num_vars == None:
        print(f"{FAIL}Number of variables is not declared correctly.{CLEAR} ({_MQ_URL})")
        return None
    system[0] = int(num_vars[1])
    num_polys = re.search(_REGEX_NUM_POLYS, contents)
    if num_polys == None:
        print(f"{FAIL}Number of polynomials is not declared correctly.{CLEAR} ({_MQ_URL})")
        return None
    system[1] = int(num_polys[1])
    if re.search(_REGEX_SYS, contents) == None:
        print(f"{FAIL}Polynomial declarations are not formatted correctly.{CLEAR} ({_MQ_URL})")
        return None
    l = 0
    for p in re.finditer(_REGEX_SYS, contents):
        system[2].append(read_poly(p[1], system[0]))
        l += 1
    if l != system[1]:
        print(f"{WARNING}Warning: File declared {system[1]} polynomials but defined {len(system[2])}{CLEAR}")
        system[1] = l
    sol = re.search(_REGEX_SOL, contents)
    if sol:
        system.append(int(sol[1]))
    else:
        print(f"{WARNING}Beware: Loaded system was not stored with a known solution.{CLEAR}")
    return system

def parse_fukuoka(file_path, print_system=False):
    print("This functionality assumes graded reverse lex order as the monomial order for all polynomials in the system.", end="\n\n")
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

def random_systems(m_low, m_high, n_low, n_high, amount):
    print(f"Generating {amount} systems...")
    systems = []
    for _ in range(amount):
        system = []
        m = randint(m_low, m_high)
        n = randint(n_low, n_high)
        ring = GF(2)[", ".join(["x" + str(i) for i in range(n)])]
        rem = 0
        for _ in range(m):
            f = ring(GF(2)[ring.gens()].random_element(degree=2, terms=Infinity))
            for v in ring.gens():
                if v^2 in f.monomials():
                    f = f + v^2 + v
            if (f != GF(2)(0)) and (f != GF(2)(1)):
                system.append(f)
            else:
                rem += 1
        m -= rem
        systems.append((system, n, m, ring, None))
    print("Finished generating systems")
    return systems

def random_systems_with_sol(m_low, m_high, n_low, n_high, amount):
    print(f"Generating {amount} systems with a known solution...")
    systems = []
    for sys_tuple in random_systems(m_low, m_high, n_low, n_high, amount):
        system, n, m, ring, _ = sys_tuple
        sol = randint(0, 2^n - 1) 
        system = [f if f(*convert(sol, n)) == 0 else (f + 1) for f in system]
        systems.append((system, n, m, ring, sol))
    print("Finished generating systems")
    return systems

def read_system(fname):
    sol = None
    sys_tuple = parse_fukuoka(fname)
    print(fname)
    n, m, system = sys_tuple[0], sys_tuple[1], sys_tuple[2]
    if len(sys_tuple) == 4:
        sol = sys_tuple[3]
    ring = GF(2)[", ".join(["x" + str(i) for i in range(n)])]
    return [(system, n, m, ring, sol)]

def fetch_systems_interactive():
    load_file = (input("Load the system from txt file? [y/n] ") == "y")
    if load_file:
        fname = input("Relative path of file to load: ")
        return read_system(fname)
    else:
        amount = int(input("Amount of systems to generate: "))
        m_tuple = tuple(int(c) for c in input("System size (specified as 'm_low m_high'): ").split(" "))
        n_tuple = tuple(int(c) for c in input("Amount of variables (specified as 'n_low n_high'): ").split(" "))
        sol_exists = (input("Ensure the system has a solution? [y/n] ") == "y")
        if sol_exists:
            return random_systems_with_sol(*m_tuple, *n_tuple, amount)
        else:
            return random_systems(*m_tuple, *n_tuple, amount)

if __name__ == "__main__":
    n, m, system = parse_fukuoka("fukuoka_test.txt", True)
    write_fukuoka("write_fukuoka_test.txt", system, n, m, GF(2)[", ".join(["x" + str(i) for i in range(n)])])
    n_1, m_1, system_1 = parse_fukuoka("write_fukuoka_test.txt", True)
    print(n == n_1, m == m_1, system == system_1)
    write_fukuoka("write_fukuokaFAIL_test.txt", system, n, m, GF(2)[", ".join(["x" + str(i) for i in range(n)])], "bad solution")

