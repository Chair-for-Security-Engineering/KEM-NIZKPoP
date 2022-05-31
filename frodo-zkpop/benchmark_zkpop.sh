#!/bin/bash

# check if all packages are installed
passflag=1
python3 -c "from natsort import natsorted" > /dev/null 2>&1
if [ $? -ne 0 ];
then
  echo "Please install Python module natsort." > /dev/stderr
  passflag=0
fi

python3 -c "from tabulate import tabulate" > /dev/null 2>&1
if [ $? -ne 0 ];
then
  echo "Please install Python module tabulate." > /dev/stderr
  passflag=0
fi

python3 -c "import glob" > /dev/null 2>&1
if [ $? -ne 0 ];
then
  echo "Please install Python module glob." > /dev/stderr
  passflag=0
fi

if [ $passflag -eq 0 ];
then
  exit
fi

echo "Frodo-640"
for i in "4 64" "8 43" "16 32" "31 26" "57 22" "107 19" "256 16" "371 15" "921 13"
do
  set -- $i
  echo "N=$1 tau=$2"
  
  make clean >> make.log 2>&1
  make OPT_LEVEL=FAST USE_OPENSSL=FALSE GENERATION_A=SHAKE128 ZKPOP_N=$1 ZKPOP_TAU=$2 NO_RESAMPLING=1 >> make.log 2>&1
  ./frodo640/test_KEM > frodo640_N$1_tau$2
done

echo "Frodo-976"
for i in "4 96" "8 64" "16 48" "31 39" "64 32" "116 28" "256 24" "424 22"
do
  set -- $i
  echo "N=$1 tau=$2"
  
  make clean >> make.log 2>&1
  make OPT_LEVEL=FAST USE_OPENSSL=FALSE GENERATION_A=SHAKE128 ZKPOP_N=$1 ZKPOP_TAU=$2 NO_RESAMPLING=1 >> make.log 2>&1
  ./frodo976/test_KEM > frodo976_N$1_tau$2
done

echo "Frodo-1344"
for i in "4 128" "8 86" "16 64" "31 52" "62 43" "107 38" "256 32" "455 29"
do
  set -- $i
  echo "N=$1 tau=$2"
  
  make clean >> make.log 2>&1
  make OPT_LEVEL=FAST USE_OPENSSL=FALSE GENERATION_A=SHAKE128 ZKPOP_N=$1 ZKPOP_TAU=$2 NO_RESAMPLING=1 >> make.log 2>&1
  ./frodo1344/test_KEM > frodo1344_N$1_tau$2
done


# create result files
python3 genfigs.py

