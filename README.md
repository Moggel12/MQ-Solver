# MQ-Solver

## TODO
[Todo-list](https://github.com/Moggel12/MQ-Solver/blob/main/TODO.md)

## Compiling code

Two targets are available via the root-directory makefile. The *default target* is `bin/mq.so` which is a shared library compiled for performance.

The `tests` target in the Makefile compiles an executable ELF file `bin/test`, compiled with sanitizers, debug symbols and other nice-to-haves. The main function, and other test-specific C-files can be found in the `test/` directory.

Compiling code using the root directory Makefile will also preparse the `src/sage/*.sage` files and store them as corresponding `src/sage/*.py` files.

While the Makefile defaults to using 32-bit integers for bitsliced polynomials and for representing the variables, this can be changed by calling make with `NSIZE=<Upper limit of bits>` and `MSIZE=<Upper limit of bits>`, with these upper limits being one of {8, 16, 32, 64, 128, 256}.

## Running tests
Both the C code and the Sage code can be run using the `run_test.py` script. This script handles running benchmarks, tests and generating polynomial systems (for testing purposes).

### Summary of `run_test.py` flags

| Short form | Long form | Purpose |
|------------|-----------|---------|
| `-l`       | `--list`  | Lists all the test functions available for invocation via this script. If a test function is not present, the function is mismatched. In order for a function to  be available, the function name must be prefixed with `test_`. Specifying this flag will cancel the usage of any other flag. |
| `-t <test name>`       | `--test <test name>`  | Runs the test function `<test name>` either using a pre-specified *MQ-challenges*-style file or generates an amount of new systems. If `<test name` is not a function from the list given by `-l` (`--list`), the script will default to running the first function listed by the list-flag. If useed with `-f` (`--file`) it will run the test function on the system from a given file, whilst not using the file-flag will give an interactive prompt for generating one or more test systems (that the test will then be run across). **If the script encounters any failed tests when running on interactively generated systems, the according failing system is dumped into a file.** |
| `-f <filename>`       | `--file <filename>`  | Specifies a file that the `-t` (`--test`) flag should load and use. Failing a test using this flag will, of course, not dump the system to a new file. This flag only has an effect when used in tandem with testing implementations |
| `-g <amount>,<n_low>,<n_high>,<m_low>,<m_high>`       | `--gen <amount>,<n_low>,<n_high>,<m_low>,<m_high>`   | Generates a directory with multiple systems, using the parameters `amount`, `n_low`, `n_high`, `m_low`, `m_high`. These parameters control the amount of systems generated as well as the lower-bounds and upper-bounds for randomly choosing `n` and `m` for each of the `amount` systems. |
| `-b <amount>,<n>,<m>` | `--bench <amount>,<n>,<m>` | Generates `amount` systems and benchmarks the `solve` function in the C-code against these. The systems are generated with `n` variables and `m` polynomials. Once the benchmarking is done, the benchmarked timings are reported to stdout. The reportings obtained are CPU times obtained through the C Time Library (`time.h`). |


### Testing code
In order to test code, whether it is the Sage code or the C code, the `run_test.py` script can be called using the `-t` (`--test`) flag (see summary of all flags above). If the flag is not specified, the script will test the first function appearing in its internal list of test functions. 

If no file is specified (i.e. no `-f` or `--file` flag specified), the script prompts the user to interactively generate systems of certain sizes. E.g. running `./run_test.py -t test_c_solve` yields the prompt (below, answers for each line was already provided):
```bash
Load the system from txt file? [y/n] n
Amount of systems to generate: 100
System size (specified as 'm_low m_high'): 10 20
Amount of variables (specified as 'n_low n_high'): 7 19
Ensure the system has a solution? [y/n] y
```
after which the script starts generating all systems and then iterates over each, running the desired test function (in the case above: `test_c_solve`). The inputs for this prompt can of course also easily specified using simple shell pipes.

If a file is specified using `-f` (`--file`), the script will run only that file against the test and no interactive prompt will appear. **Currently, the script only supports loading one test-file at a time**.

The functions available for testing can be found by using the `-l` (`--list`) flag.

Whenever a test fails on an interactively generated file, the necessary information is dumped into an *MQ-challenges*-like file that can be easily loaded using `-f` (`--file`) in later calls. The format of the dumped file is `<digest of file>_<name of failed test function>_system.txt`.

**Running tests using `-t` (`--test`), or running the script with no flags, requires that the default make target was compiled and that a `bin/mq.so` is present in the project root.**

### Benchmarking code

Using the `-b` (`--bench`) flag, the C code can be benchmarked end-to-end (i.e. benchmarking the solver from initiation to successfully finding a solution). **Currently, benchmarking the Sage code is not supported via `run_test.py` script, and can only be done manually in the appropriate `.sage` files**.

The `-b` (`--bench`) flag requires three comma-separated inputs `<amount>,<n>,<m>` that tells the script how many systems the benchmarking should be done through (timings are averaged in the end) as well as the size of the system (`n` being the amount of variables and `m` being the amount of polynomials).

### Generating tests

Test cases can be generated using the `-g` (`--gen`) flag, which takes as arguments five comma-separated values: `amount`, `n_low`, `n_high`, `m_low` and `m_high`. These values tell the script how many systems to generate as well as the lower- and upper-bounds for randomly assigning amount of polynomials and variables.

The systems generated this way reside in the `test_systems` directory. The format for the names of these files is `system_<n>_<m>_<digest of file>.txt`.
