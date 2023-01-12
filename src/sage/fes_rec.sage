from itertools import combinations
from fes import update, fes_eval
import numpy as np
from utils import convert, bitslice
import time

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

def fes_recover(system, n, n1, degree, ring):    # parameter f for debugging polynomial interpolation
    res = [None] * 2^(n - n1)
    s = None
    prefix = []
    d = {}

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


    print(fes_time_eval, fes_time_inter)

    return res

def test():
    n = 8
    degree = 3

    ring = GF(2)[", ".join(["x" + str(i) for i in range(n)])]

    print(ring)
    print()


    X = R.gens()

    f = GF(2).random_element()

    for d in range(1, degree+1):
      for mon in combinations(range(n), d):
        f += prod([X[i] for i in mon]) * GF(2).random_element()
      

    print("Input polynomial:", f)
    print()

    print("Performing full evaluation...")
    eval_full = fes_eval(f)

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
    test()
