#! /usr/bin/env python3

import argparse
import time
import sys
import os
import hashlib

from utils import fetch_systems_interactive, read_system, SUCCESS, CLEAR, FAIL, bitslice, random_systems_with_sol
from c_config import TEST_BIN_AVAILABLE
from fes_rec import *
from dinur import *
from fes import *
from mob_new import *

_test_functions = {key: val for key, val in globals().items() if callable(val) and key.startswith("test_")}

_availability_filter = lambda f_name : ("_SAN" in f_name) if TEST_BIN_AVAILABLE else ("_SAN" not in f_name)

_available_tests = {key: val for key, val in _test_functions.items() if _availability_filter(key)}
_unavailable_tests = {key: val for key, val in _test_functions.items() if key not in _available_tests}
_default_test = list(_test_functions.items())[0]

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run certain tests via the commandline (C or Sage/Python). Beware, this tool will not sanitize any dirty inputs.")
    parser.add_argument("-f", "--file", help="Load an MQ-challenges style file and use it for testing purposes", type=str)
    parser.add_argument("-t", "--test", help="Run a test based on the names listed via --list (-l)", type=str)
    parser.add_argument("-l", "--list", help="List all available test functions", default=False, action="store_true")
    parser.add_argument("-g", "--gen", help="Generate MQ-challenge style test files. Expects a comma-separated input: '<amount>,<n_low>,<n_high>,<m_low>,<m_high>'", type=str)
    parser.add_argument("-b", "--bench", help="Measure how much the C code can bench. Expects a comma-separated input: '<amount>,<n>,<m>'", type=str)
    return parser.parse_args()

def list_tests():
    print(f"{SUCCESS}Tests that may be run:{CLEAR}")
    for func_name, _ in _available_tests.items():
        print(f"\t{SUCCESS}{func_name}{CLEAR}")
    if TEST_BIN_AVAILABLE:
        print(f"{FAIL}The following tests are unavailable as they require the bin/mq.so file to run (see README for compilation instructions){CLEAR}")
    else:
        print(f"{FAIL}The following tests are unavailable as they require the bin/test file to run (see README for compilation instructions){CLEAR}")
    for func_name, _ in _unavailable_tests.items():
        print(f"\t{FAIL}{func_name}{CLEAR}")

def call_test(all_tuples, test, write_file):
    print(all_tuples)
    succeeded = False
    for i, sys_tuple in enumerate(all_tuples):
        print(f"\n{ITER}== {i + 1} =={CLEAR}")
        print(sys_tuple)
        succeeded = test[1](sys_tuple)
        if not succeeded and write_file:
            m = hashlib.sha256()
            m.update(bytes(str(sys_tuple), "utf-8"))
            write_fukuoka(f"{m.hexdigest()}_failed_{test[0]}_system.txt", *sys_tuple)
            break
    if succeeded: 
        print(f"{SUCCESS}Found no errors in {len(all_tuples)} trials{CLEAR}", file=sys.stderr)
    else:
        print(f"{FAIL}Error found in <{test[0]}>{CLEAR}", file=sys.stderr)

def gen_files(amount, n_low, n_high, m_low, m_high):
    gen_dir = os.path.join(os.getcwd(), "test_systems")
    if not os.path.exists(gen_dir):
        os.mkdir(gen_dir)
    for sys_tuple in random_systems_with_sol(m_low, m_high, n_low, n_high, amount):
        m = hashlib.sha256()
        m.update(bytes(str(sys_tuple), "utf-8"))
        write_fukuoka(f"{gen_dir}/system_{sys_tuple[1]}_{sys_tuple[2]}_{m.hexdigest()}.txt", *sys_tuple)

def main():
    all_tuples = None
    test = _default_test
    write_file = True
    args = parse_arguments()
    if args.list:
        list_tests()
        return
    if args.gen:
        gen_files(*[int(str_num) for str_num in args.gen.split(",")[:5]])
        return
    elif args.bench:
        amnt, n, m = [int(str_num) for str_num in args.bench.split(",")[:3]]
        system_tuples = random_systems_with_sol(m, m, n, n, amnt)
        c_benchmark(system_tuples)
        return
    if not args.file:
        all_tuples = fetch_systems_interactive()
    else:
        all_tuples = read_system(args.file)
        write_file = False
    if args.test:
        if args.test in _test_functions:
            test = (args.test, _test_functions[args.test])
        else:
            print(f"{WARNING}Could not recognize test function. Using default test ({_default_test}){CLEAR}", file=sys.stderr)
    call_test(all_tuples, test, write_file)


if __name__ == "__main__":
    main()