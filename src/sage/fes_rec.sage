from itertools import combinations
from fes import update, fes_eval, c_bruteforce
import numpy as np
from utils import convert, bitslice, fetch_c_func
import time

import ctypes as ct
from c_config import *

from utils import random_systems_with_sol

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
    res[0] = new_parities 
    d[0] = new_parities 

    fes_time_eval = 0
    fes_time_inter = 0

    for si in range(1, 2^(n - n1)):
        if len(bits(si)) > degree:
            # We have the required derivatives; compute the missing evluation value.
            fes_time_eval -= time.time()

            k = bits(si)[:degree]

            for j in reversed(range(0, len(k))):
                d[sum([2^i for i in k[:j]])] = int(d[sum([2^i for i in k[:j]])]) ^^ int(d[sum([2^i for i in k[:j+1]])])
        

            fes_time_eval += time.time()

        else:
            # We need to interpolate derivatives.
            fes_time_inter -= time.time()

            k = bits(si)[:degree]
    
            prefix = [pos for pos,b in enumerate(reversed(bin(si ^^ (si >> 1))[2:])) if b == "1"]

            s, new_parities = part_eval(system, prefix, n, n1, s)

            fes_time_inter += time.time()

            fes_time_eval -= time.time()

            prev = d[0]
            d[0] = new_parities

            for j in range(1, len(k)+1):

                if j < len(k):
                    tmp = d[sum([2^i for i in k[:j]])]

                d[sum([2^i for i in k[:j]])] = int(d[sum([2^i for i in k[:j-1]])]) ^^ int(prev)
    
                if j < len(k):
                    prev = tmp

            fes_time_eval += time.time()
        
        res[si ^^ (si >> 1)] = d[0]

    return res

def c_fes_recover(system, n, n1, deg):
    c_results = (C_VARS_T * int(2^(n - n1)))()
    c_system = (C_POLY_T * len(system))(*system)

    args = [Type.P(C_POLY_T), Type.U, Type.U, Type.U, Type.P(C_VARS_T)]
    res = Type.U8
    recover = fetch_c_func("fes_recover", args, res)
    error = recover(c_system, n, n1, deg, c_results)

    if error == 1:
        return None

    py_list = []
    for i in range(2^(n - n1)):
        py_list.append(int(c_results[i]))

    return py_list

def test_c_fes_recover(sys_tuple):
    system, n, _, ring, _ = sys_tuple
    n1 = int(ceil(n/(5.4))) # Quadratic systems are assumed here, see page 19 of full dinur paper for explanation
    d = randint(1, n - n1)
    print(system)
    print(n, n1, d)
    system = bitslice(system, ring.gens())
    c_results = c_fes_recover(system, n, n1, d + 1)
    if c_results == None:
        print("Error with memory allocation in C code.")
        return False
    py_results = fes_recover(system, n, n1, d + 1)
    if len(c_results) != len(py_results):
        print("Solutions are not of equal length")
        print(c_results)
        print(py_results)
        return False
    for c_r, p_r in zip(c_results, py_results):
        if (c_r != p_r):
            print("Solutions differ!", c_r, "!=", p_r)
            print(c_results, py_results)
            return False
    return True
if __name__ == "__main__":
    # test()
    # fes_recover(sl_sys, n, n1, d)
    test_c_fes_recover(1, 7, 7)
