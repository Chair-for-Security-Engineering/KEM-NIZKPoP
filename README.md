# Proof-of-possession for KEM certificates using verifiable generation

This archive contains the artifact for the paper "Proof-of-possession for KEM certificates using verifiable generation".

Specifically, it contains:

- Python and Sage scripts for estimating proof size parameters based on formulas in the paper.
- a C implementation of the verifiable generation procedure for FrodoKEM
- a C implementation of the verifiable generation procedure for Kyber
- benchmarking and profiling scripts

The implementation is based on the original C implementations of [FrodoKEM](https://github.com/microsoft/PQCrypto-LWEKE) and [Kyber](https://github.com/pq-crystals/kyber).

## Source code organization

- `frodo-zkpop`: The FrodoKEM implementation with verifiable generation, based on the [original FrodoKEM implementation](https://github.com/microsoft/PQCrypto-LWEKE). Files of interest specifically for the verifiable generation modifications include:
	- `frodo-zkpop/src/zkpop.c`: The main code for our zero-knowledge verifiable generation routine.
	- `frodo-zkpop/tests/test_kem.c`: Test harness and benchmarking script adapted to include our proof generation and verification.
	- `frodo-zkpop/benchmark_zkpop.sh`: Script to run the verifiable generation benchmarks.
- `kyber-zkpop`: The Kyber implementation with verifiable generation, based on the [original Kyber implementation](https://github.com/pq-crystals/kyber). Files of interest specifically for the verifiable generation modifications include:
	- `kyber-zkpop/avx2/zkpop.{h,c}`: The main code for our zero-knowledge verifiable generation routine.
	- `kyber-zkpop/avx2/test_zkpop.c`: Test harness for our proof generation and verification.
	- `kyber-zkpop/avx2/test_speed.c`: Benchmarking script adapted to include our proof generation and verification.
	- `kyber-zkpop/avx2/benchmark_zkpop.sh`: Script to run the verifiable generation benchmarks.
- `benchmark.sh`: Top-level script to run all the verifiable generation benchmarks.
- `scripts`: Scripts to calculate proof parameters and ideal proof sizes.
	- `scripts/calculate_M_{frodo,kyber}.sage`: Sage script to calculate gamma and M values.
	- `scripts/params.py`: Python3 script to calculate ideal proof sizes for different parameter regimes.

## Requirements

Both implementations require a processor with support for AVX2.

The software has been built and tested on Arch Linux using gcc 11.2.0, running on an Intel processor with AVX2 extensions.

To compile, ensure you satisfy the build requirements for FrodoKEM and Kyber, described in their respective README files.  At the minimum, you need a recent version of gcc or clang.

In order to do operation profiling, you must have gprof installed.

Python3 is required.  The following Python packages must be installed:

```bash
# Using apt on Debian/Ubuntu:
apt install python3-natsort python3-tabulate
# Using pip3:
pip3 install natsort scipy tabulate
```

Sage is required to run the `calculate_M_*.sage` scripts.

## Instructions

**To obtain performance data:**

- In the root directory of the folder, run `./benchmark.sh` to obtain performance data.
- In `kyber-zkpop/avx2`, run `make clean && make test_speed512 ZKPOP_N=65536 ZKPOP_TAU=8 NO_RESAMPLING=1 && ./test_speed512` to obtain the benchmark for N=65536, tau=8
- In `frodo-zkpop`, run `make clean && make OPT_LEVEL=FAST USE_OPENSSL=FALSE GENERATION_A=SHAKE128 ZKPOP_N=65536 ZKPOP_TAU=8 && frodo640/test_KEM` to obtain the benchmark for N=65536, tau=8

This will run the Kyber and FrodoKEM verifiable generation benchmarks at each of three security levels, and at several (N, tau) parameter sizes, with 100 iterations in each configuration.

- Raw FrodoKEM results will be saved in the files `frodo-zkpop/frodo{640,976,1344}_N*_tau*`
- A table summarizing FrodoKEM results will be saved in the files `frodo-zkpop/frodocycles{512,768,1024}`
- Raw Kyber results will be saved in the files `kyber-zkpop/avx2/kyber{512,768,1024}_N*_tau*`
- A table summarizing Kyber results will be saved in the files `kyber-zkpop/avx2/kybercycles{512,768,1024}`

If any file is empty, please check `kyber-zkpop/avx2/make.log` or `frodo-zkpop/make.log`.

**To obtain profiling data:**

Ensure that `gprof` is available.

- In `kyber-zkpop/avx2`, run `make clean && make test_profile ZKPOP_N=4 ZKPOP_TAU=64 NO_RESAMPLING=1 && ./test_profile`, then `gprof ./test_profile` to obtain the profiling data.

**To generate proof parameters and ideal proof sizes:**

- To calculate gamma and M for FrodoKEM and Kyber:
	- In the `scripts` directory, run `sage calculate_M_frodo.py` and `sage calculate_M_kyber.py`.
- To calculate ideal proof sizes for the various parameter regimes:
	- In the `scripts` directory, run `python3 params.py`.
