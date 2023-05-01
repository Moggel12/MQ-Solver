import ctypes as ct
import sys
import argparse
import os
import math

# from utils import parse_fukuoka, write_fukuoka, random_systems, random_systems_with_sol, fetch_c_func

_libc = ct.CDLL("libc.so.6")

RSEED = 42
MAX_HISTORY = 30

TEST_BIN_AVAILABLE = os.path.exists(os.path.join(os.path.dirname(__file__), "../../bin/test"))

C_POLY_T = None
C_VARS_T = None

C_VECTORIZED = False
C_VECTOR_SIZE = 0
C_FIXED_VARS = 0

# Typedefs
class Type():
    U   = ct.c_uint
    U8  = ct.c_uint8 
    U16 = ct.c_uint16
    U32 = ct.c_uint32
    U64 = ct.c_uint64 
    SZ = ct.c_size_t
    P8  = ct.POINTER(ct.c_uint8)
    P16 = ct.POINTER(ct.c_uint16)
    P32 = ct.POINTER(ct.c_uint32)
    P64 = ct.POINTER(ct.c_uint64)

    def P(c_type):
        return ct.POINTER(c_type)

type_dict = {
    8:   Type.U8,
    16:  Type.U16,
    32:  Type.U32,
    64:  Type.U64,
    128: Type.U64,
    256: Type.U64
}
try:
    with open(os.path.join(os.path.dirname(__file__), ".compile_config"), "r") as f:
        bits_str = f.readline().split(" ")
        bits = int(bits_str[0])
        if bits >= 128:
            C_VECTORIZED = True
            C_VECTOR_SIZE = bits//int(bits_str[1])
            C_FIXED_VARS = int(math.log2(C_VECTOR_SIZE)) - 2
            bits = 64
        C_POLY_T = type_dict[bits]
        C_VARS_T = type_dict[bits]

except FileNotFoundError:
    C_POLY_T = Type.U32
    C_VARS_T = Type.U32

def srand(seed):
    _libc.srand(seed)

def rand():
    return _libc.rand()


