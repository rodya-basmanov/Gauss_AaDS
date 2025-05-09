# Vectorized Gauss Elimination Method

Implementation of the Gauss elimination method for solving systems of linear algebraic equations using the Eigen library for vectorized operations.

## Requirements

- C++17 compatible compiler
- Eigen library
- Google Test (for running tests)
- OpenMP (for parallel computing)

## Building and Running

To build the program:

```
make
```

Or using the script:

```
chmod +x build.sh
./build.sh
```

Running tests:

```
chmod +x tests.sh
./tests.sh
```

## Usage

### Reading system from CSV file

```
./main input.csv output.csv
```

### Generating random system

```
./main --generate <size> <seed> <output_file>
```

Examples:

```
./main example.csv solution.csv
./main --generate 10 42 random_system.csv
```

## CSV File Format

### Input File

Coefficient matrix A and right-hand side vector b in the form [A|b]:

```
a11,a12,...,a1n,b1
a21,a22,...,a2n,b2
...
am1,am2,...,amn,bm
```

### Output File

Solution in the form:

```
x0,<value>
x1,<value>
...
xn-1,<value>
```