import ctypes as ct
from itertools import combinations
from fes import update, fes_eval
import numpy as np
from utils import convert, bitslice, fetch_c_func
import time

from utils import random_system_with_sol

# Get index position of first bit set (if any).
def bit1(x):
    x = x&-x

    x = int(x).bit_length()

    return None if x == 0 else x-1

# Get indices of set bits as array.
def bits(x):
    if x == 0:
     return None

    ret = []
    while x > 0: 
        ret.append(bit1(x))

        x = x ^^ (x&-x)

    return ret

def part_eval(system, prefix, n, n1, s):
    s = update(s, system, n, n1, prefix)
    U_parities = fes_eval(system, n, n1, prefix, s, True)
    return s, U_parities

def fes_recover(system, n, n1, degree, ring=None):    # parameter f for debugging polynomial interpolation
    res = [None] * 2^(n - n1)
    s = None
    prefix = []
    d = {}

    if ring != None:
        system = bitslice(system, ring.gens())

    s, new_parities = part_eval(system, prefix, n, n1, s)
    print("Got parities", new_parities)
    res[0] = new_parities 
    d[0] = new_parities 

    fes_time_eval = 0
    fes_time_inter = 0

    for si in range(1, 2^(n - n1)):
        if len(bits(si)) > degree:
            # print("=> Compute evaluation")
            # We have the required derivatives; compute the missing evluation value.
            fes_time_eval -= time.time()

            k = bits(si)[:degree]

            for j in reversed(range(0, len(k))):
                # print(f"For si = {si} = {bits(si)[:degree]} we have indices: {sum([2^i for i in k[:j]])}, {sum([2^i for i in k[:j+1]])}")
                d[sum([2^i for i in k[:j]])] = int(d[sum([2^i for i in k[:j]])]) ^^ int(d[sum([2^i for i in k[:j+1]])])
        

            fes_time_eval += time.time()

        else:
            # print("=> Interpolate derivative")
            # We need to interpolate derivatives.
            fes_time_inter -= time.time()

            k = bits(si)[:degree]
    
            prefix = [pos for pos,b in enumerate(reversed(bin(si ^^ (si >> 1))[2:])) if b == "1"]

            s, new_parities = part_eval(system, prefix, n, n1, s)
            print("Got parities", new_parities)

            fes_time_inter += time.time()

            fes_time_eval -= time.time()

            prev = d[0]
            d[0] = new_parities

            for j in range(1, len(k)+1):
                # print(f"Computed so far: {[0]}")

                if j < len(k):
                    tmp = d[sum([2^i for i in k[:j]])]
    
                # print(f"For si = {si} = {bits(si)[:degree]} we have indices: {sum([2^i for i in k[:j]])}, {sum([2^i for i in k[:j-1]])}")
                d[sum([2^i for i in k[:j]])] = int(d[sum([2^i for i in k[:j-1]])]) ^^ int(prev)
    
                if j < len(k):
                    prev = tmp

            fes_time_eval += time.time()
        
        res[si ^^ (si >> 1)] = d[0]


    # print(fes_time_eval, fes_time_inter)

    return res

def c_fes_recover(system, n, n1, deg):
    c_results = (ct.c_uint8 * int(2^(n - n1)))()
    c_system = (ct.c_uint8 * len(system))(*system)

    args = [ct.POINTER(ct.c_uint8), ct.c_uint, ct.c_uint, ct.c_uint, ct.POINTER(ct.c_uint8)]
    res = ct.c_uint8
    recover = fetch_c_func("fes_recover", args, res)
    recovered = recover(c_system, n, n1, deg, c_results)

    if recovered == 1:
        return None

    py_list = []
    for i in range(2^(n - n1)):
        py_list.append(int(c_results[i]))

    return py_list
        

def test_c_fes_recover(trials, m = 5, n = 5):
    for _ in range(1):
        ring.<x0,x1,x2,x3,x4> = GF(2)[]
        system = [x0*x2 + x2*x3 + x0*x4 + x2*x4 + x0 + x2 + x3, x0*x2 + x1*x2 + x0*x3 + x1*x3 + x2*x3 + x0*x4 + x1*x4 + x2*x4 + x1, x0*x1 + x0*x2 + x1*x2 + x0*x3 + x1*x3 + x2*x3 + x0*x4 + x1*x4 + x3*x4 + x0 + x1 + x2 + x3 + x4 + 1, x0*x1 + x0*x2 + x1*x3 + x2*x3 + x1*x4 + x3*x4 + x0 + x1 + x2, x0*x2 + x0*x4 + x1*x4 + x3*x4 + x2 + x3]
        n = 5
        # system, _, ring, n = random_system_with_sol(m, m, n, n)
        n1 = int(ceil(n/(5.4))) # Quadratic systems are assumed here, see page 19 of full dinur paper for explanation
        # d = randint(1, n - n1)
        d = 3
        print(system)
        print(n, n1, d)
        print("\n### C ###\n")
        system = bitslice(system, ring.gens())
        c_results = c_fes_recover(system, n, n1, d + 1)
        if c_results == None:
            print("Error with memory allocation in C code.")
            return 
        print("\n### Python ###\n")
        py_results = list(fes_recover(system, n, n1, d + 1, ring))
        if len(c_results) != len(py_results):
            print("Solutions are not of equal length")
            print(c_results)
            print(py_results)
            return
        for c_r, p_r in zip(c_results, py_results):
            if (c_r != p_r):
                print("Solutions differ!", c_r, "!=", p_r)
                print(c_results, py_results)
                return
    print("Found no errors")

def test():
    n = 8
    degree = 3

    ring = GF(2)[", ".join(["x" + str(i) for i in range(n)])]

    print(ring)
    print()


    X = ring.gens()

    f = GF(2).random_element()

    for d in range(1, degree+1):
      for mon in combinations(range(n), d):
        f += prod([X[i] for i in mon]) * GF(2).random_element()
      

    print("Input polynomial:", f)
    print()

    print("Performing full evaluation...")
    eval_full = fes_eval(f, n)

    print(eval_full)

    check = [f(*reversed([int(i) for i in f"{x:0{n}b}"])) for x in range(2^n)]

    print(check)

    print(eval_full == check)


    print("\nPerforming partial evaluation and recovery...")
    eval_part = fes_eval(f, f.degree())

    print("[" + ", ".join(['-' if v == None else str(v) for v in eval_part]) + "]")

    recovered = fes_recover(eval_part, len(f.parent().gens()), f.degree(), f, ring)

    print(recovered)

    print(recovered == eval_full)

if __name__ == "__main__":
    # test()
    test_c_fes_recover(20)
