# MQ-Solver

## Compiling the report

The report PDF can be compiled by running `make pdf`, and can be cleaned up by running `make pdfclean`. 

Note: The LaTeX code for the report uses the `minted` package, which expects `pygments` package (or `python-pygments`) to be installed on the system.

## Compiling code

Two targets are available via the root-directory makefile. The *default target* is `bin/mq.so` which is a shared library compiled for performance.

The `tests` target in the Makefile compiles an executable ELF file `bin/test`, compiled with sanitizers, debug symbols and other nice-to-haves. The main function and other test-specific C-files can be found in the `test/` directory.

Compiling code using the root directory Makefile will also preparse the `src/sage/*.sage` files and store them as corresponding `src/sage/*.py` files.

While the Makefile defaults to using 32-bit integers for bit-sliced polynomials and for representing the variables, this can be changed by calling make with `REGSIZE=<Upper limit of bits>`, with these upper limits being one of {8, 16, 32, 64, 128, 256}. Specifying 128 or 256 bits implies that AVX or AVX2 instructions, respectively, are available on the CPU.

## Docker

The repository also contains a Dockerfile with an Ubuntu image and all the necessary packages installed to compile and run this project. Simply build an image from this Dockerfile and share the repository with the container!  

## Executing the code
Both the C code and the Sage code can be run using the `run.py` script. This script handles running benchmarks, tests, and generating polynomial systems (for testing purposes). The script is also a front-end to the solver. If the *shared library* (`bin/mq.so`) has been compiled, the script can run a `-O3` optimized solver, essentially acting as a front-end to the shared library.

To use the script as a front-end to the solver(s), simply run `./run.py` without arguments to get an interactive prompt for generating input systems, or use the `-f` (`--file`) flag to point the solver to an [MQ-challenge](https://www.mqchallenge.org/)-styled file used as input. The flags `-F` (`--FES`) and `-p` (`--parallelize`) will alter the behavior of the solver before trying to solve the system. See the table below for more information on these flags.

### Summary of `run.py` flags

| Short form | Long form | Purpose |
|------------|-----------|---------|
| `-l`       | `--list`  | Lists all the test functions available for invocation via this script. If a test function is not present, the function is mismatched. In order for a function to  be available, the function name must be prefixed with `test_`. Specifying this flag will cancel the usage of any other flag. |
| `-t <test name>`       | `--test <test name>`  | Runs the test function `<test name>` either using a pre-specified *MQ-challenges*-style file or generates an amount of new systems. If `<test name` is not a function from the list given by `-l` (`--list`), the script will default to running the first function listed by the list-flag. If useed with `-f` (`--file`) it will run the test function on the system from a given file, whilst not using the file-flag will give an interactive prompt for generating one or more test systems (that the test will then be run across). **If the script encounters any failed tests when running on interactively generated systems, the according failing system is dumped into a file.** |
| `-f <filename>`       | `--file <filename>`  | Specifies a file that the `-t` (`--test`) flag should load and use. Failing a test using this flag will, of course, not dump the system to a new file. This flag only has an effect when used in tandem with testing implementations |
| `-g <amount>,<n_low>,<n_high>,<m_low>,<m_high>`       | `--gen <amount>,<n_low>,<n_high>,<m_low>,<m_high>`   | Generates a directory with multiple systems, using the parameters `amount`, `n_low`, `n_high`, `m_low`, `m_high`. These parameters control the amount of systems generated as well as the lower-bounds and upper-bounds for randomly choosing `n` and `m` for each of the `amount` systems. |
| `-b <amount>,<n>,<m>` | `--bench <amount>,<n>,<m>` | Generates `amount` systems and benchmarks the `solve` function in the C-code against these. The systems are generated with `n` variables and `m` polynomials. Once the benchmarking is done, the benchmarked timings are reported to stdout. The reportings obtained are CPU times obtained through the C Time Library (`time.h`). |
| `-F` | `--FES` | Run the built-in (in C) FES procedure either using the `-b` (`--bench`) flag, or purely as a solver (parallelized *and* sequential are available). The presence of this flag switches these actions to use FES instead of Dinur's solver. This flag does not work with the `-t` (`--test`) flag. |
| `-p` | `--parallelize` | Can be used when running the script as an MQ-solver on some input system. This flag tells the script to spawn multiple independent sub-processes to solve the input system, by fixing variables and dividing the new fixed systems up between the sub-processes. Using this flag will spawn a power-of-two amount of sub-processes, specifically depending on the available cores on the host. Currently, the script only terminates when *all* spawned processes have finished, however, it will print the solution found when a sub-process terminates. |


### Testing code
In order to test code, whether it is the Sage code or the C code, the `run.py` script can be called using the `-t` (`--test`) flag (see the summary of all flags above). If the flag is not specified, the script will test the first function appearing in its internal list of test functions. 

If no file is specified (i.e. no `-f` or `--file` flag specified), the script prompts the user to interactively generate systems of certain sizes. E.g. running `./run.py -t test_c_solve` yields the prompt (below, answers for each line were already provided):
```bash
Load the system from txt file? [y/n] n
Amount of systems to generate: 100
System size (specified as 'm_low m_high'): 10 20
Amount of variables (specified as 'n_low n_high'): 7 19
Ensure the system has a solution? [y/n] y
```
after which the script starts generating all systems and then iterates over each, running the desired test function (in the case above `test_c_solve`). The inputs for this prompt can of course also be easily specified using simple shell pipes.

If a file is specified using `-f` (`--file`), the script will run only that file against the test and no interactive prompt will appear. **Currently, the script only supports **loading one **test file** at a** time**.

The functions available for testing can be found by using the `-l` (`--list`) flag.

Whenever a test fails on an interactively generated file, the necessary information is dumped into an *MQ-challenges*-like file that can be easily loaded using `-f` (`--file`) in later calls. The format of the dumped file is `<digest of file>_<name of failed test function>_system.txt`.

**Running some tests using `-t` (`--test`), or running the script as solver front-end, requires that the default make target was compiled (potentially with a different `REGSIZE`) and that a `bin/mq.so` is present in the project root.**

### Benchmarking code

Using the `-b` (`--bench`) flag, the C code can be benchmarked end-to-end (i.e. benchmarking the solver from initiation to successfully finding a solution). **Currently, benchmarking the Sage code is not supported via `run.py` script, and can only be done manually in the appropriate `.sage` files**.

The `-b` (`--bench`) flag requires three comma-separated inputs `<amount>,<n>,<m>` that tell the script how many systems the benchmarking should be done through (timings are averaged in the end) as well as the size of the system (`n` is the number of variables and `m` is the number of polynomials). Running this will generate a directory `bench_<n>_<m>_<time>/` containing the systems used for the benchmark and a `bench.csv` file with most of the benchmark values present. The benchmark values that are dumped to the CSV file are specified *by name* in the `C_BENCHMARK_VARS` global, in `src/sage/c_config.py`.

**Alternatively**, the benchmarking procedure can be pointed to a directory generated by the approach detailed above by specifying `-b <directory>` instead of `-b <amount>,<n>,<m>`. This will automatically update the CSV file in the directory (or create it) if necessary. Be sure that the columns in the CSV file stay the same between two benchmarks using the same directory!

### Generating tests

Test cases can be generated using the `-g` (`--gen`) flag, which takes as arguments five comma-separated values: `amount`, `n_low`, `n_high`, `m_low` and `m_high`. These values tell the script how many systems to generate as well as the lower- and upper-bounds for randomly assigning amounts of polynomials and variables.

The systems generated this way reside in the `test_systems` directory. The format for the names of these files is `system_<n>_<m>_<digest of file>.txt`.
