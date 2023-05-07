#! /usr/bin/env python3

import argparse
import time
import sys
import os
import hashlib
import multiprocessing
import math as m

from src.sage.utils import fetch_systems_interactive, read_system, SUCCESS, CLEAR, FAIL, bitslice, random_systems_with_sol, fix_variables
from src.sage.c_config import TEST_BIN_AVAILABLE
from src.sage.fes_rec import *
from src.sage.dinur import *
from src.sage.fes import *
from src.sage.mob_new import *

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
    parser.add_argument("-p", "--parallelize", help="Spawns parallel processes according to the available threads on the CPU, each with a set of variables fixed. This flag only affects using the script for solving systems.", default=False, action="store_true")
    return parser.parse_args()

def solve(all_tuples):
    print(f"Solving {len(all_tuples)} systems one-by-one..")
    if TEST_BIN_AVAILABLE:
        print(f"{FAIL}Compile shared object with the Makefile to use this flag.{CLEAR}")
        return None
    for i, sys_tuple in enumerate(all_tuples):
        print(f"\n{ITER}== {i + 1} =={CLEAR}")
        system, n, m, ring, _ = sys_tuple
        sl_system = bitslice(system, ring.gens())
        sol = c_solve(sl_system, n, m)
        if not sol:
            print("=> C code did not solve the system")
        else:
            print(f"=> Found solution: {sol}")

def p_solve(challenge):
    sys_tuple, fixture = challenge
    system, n, m, ring, _ = sys_tuple
    sl_system = bitslice(system, ring.gens())
    start_time = time.perf_counter_ns()
    solution = c_solve(sl_system, n, m)
    elapsed_time = time.perf_counter_ns() - start_time
    print(f"Process found {solution + fixture if solution else None} in {round(elapsed_time / 1000000)}ms")

def parallelize(all_tuples):
    core_count = multiprocessing.cpu_count()
    print(f"Running {core_count} instances in parallel..")
    if TEST_BIN_AVAILABLE:
        print(f"{FAIL}Compile shared object with the Makefile to use this flag.{CLEAR}")
        return None
    fixed_vars = int(m.log2(core_count))
    challenges = []
    with multiprocessing.Pool((1**fixed_vars)) as pool:
        for sys_tuple in all_tuples:
            challenges = [(fixed_tuple, convert(i, fixed_vars)) for i, fixed_tuple in enumerate(fix_variables(sys_tuple, fixed_vars))]
            pool.map(p_solve, challenges)

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
    print(f"Running {len(all_tuples)} tests")
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
        if args.parallelize:
            print(f"{WARNING}Tests will not run parallelized, however, the -p (--p) flag was given.{CLEAR}")
        test = _default_test
        if args.test in _test_functions:
            test = (args.test, _test_functions[args.test])
        else:
            print(f"{WARNING}Could not recognize test function. Using default test ({_default_test}){CLEAR}", file=sys.stderr)
        call_test(all_tuples, test, write_file)
        return
    if args.parallelize:
        parallelize(all_tuples)
    else:
        solve(all_tuples)


if __name__ == "__main__":
    main()
