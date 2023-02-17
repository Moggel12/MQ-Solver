#! /usr/bin/env python3

import argparse
import time

from utils import fetch_systems_interactive, read_system, SUCCESS
from fes_rec import *
from dinur import *
from fes import *
from mob_new import *

_test_functions = {key: val for key, val in globals().items() if callable(val) and key.startswith("test_")}
_default_test = list(_test_functions.items())[0]

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run certain tests via the commandline (C or Sage/Python).")
    subparser = parser.add_subparsers()
    parser.add_argument("-f", "--file", help="Load an MQ-challenges style file and use it for testing purposes", type=str)
    parser.add_argument("-t", "--test", help="Run a test based on the names listed via --list (-l)", type=str)
    parser.add_argument("-l", "--list", help="List all available test functions", default=False, action="store_true")
    return parser.parse_args()

def list_tests():
    for func_name, _ in _test_functions.items():
        print(func_name)

def call_test(all_tuples, test, write_file):
    succeeded = False
    for sys_tuple in all_tuples:
        succeeded = test[1](sys_tuple)
        if not succeeded and write_file:
            write_fukuoka(f"{int(time.time())}_failed_{test[0]}_system.txt", *sys_tuple)
            break
    if succeeded: 
        print(f"{SUCCESS}Found no errors in {len(all_tuples)} trials{CLEAR}")
    else:
        print(f"{FAIL}Error found in <{test[0]}>{CLEAR}")

def main():
    all_tuples = None
    test = _default_test
    write_file = True
    args = parse_arguments()
    if args.list:
        list_tests()
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
            print(f"{WARNING}Could not recognize test function. Using default test ({_default_test}){CLEAR}")
    call_test(all_tuples, test, write_file)


if __name__ == "__main__":
    main()