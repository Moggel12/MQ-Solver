import ctypes as ct
import sys
import argparse

_SRC_DIR = "../src/sage"

sys.path.insert(0, _SRC_DIR)

from utils import parse_fukuoka, write_fukuoka, random_systems, random_systems_with_sol, fetch_c_func

_libc = ct.CDLL("libc.so.6")

RSEED = 0

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

C_POLY_T = Type.U16
C_VARS_T = Type.U16

def srand(seed):
    _libc.srand(seed)

def rand():
    return _libc.rand()

def call_c_test(name, args, args_type, res_type):
    return fetch_c_func(name, args_type, res_type, True)(*args)


