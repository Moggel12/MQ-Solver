#! /usr/bin/env python3

import sys
from math import comb

if len(sys.argv) != 3:
    print("Please provide the dimensions of the lookup table to be generated.")
    sys.exit(0)

if (not sys.argv[1].isdigit()) or (not sys.argv[2].isdigit()):
    print("Arguments must be integers.")
    sys.exit(0)

n = int(sys.argv[1])
m = int(sys.argv[2])

prefix_code = (
    f"#ifndef BINOM_H\n"
    f"#define BINOM_H\n"
     "\n"
    f"#define BINOM_DIM1 {n}\n"
    f"#define BINOM_DIM2 {m}\n"
     "\n"
    f"unsigned int lk_binom[{n * m}] = \n"
     "{"
)

postfix_code = (
    "};\n"
    "\n"
    "#endif  // !BINOM_H"
)

print(prefix_code)
for i in range(n):
    print(end=" "*2)
    for j in range(m):
        print(f"{comb(i, j)}u", end=", ")
    print()
print(postfix_code)