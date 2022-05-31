#!/bin/bash

cd kyber-zkpop/avx2
bash benchmark_zkpop.sh

cd ../../frodo-zkpop
bash benchmark_zkpop.sh
