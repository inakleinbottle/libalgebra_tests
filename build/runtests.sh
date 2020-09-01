#! /bin/bash
set -e

FNAME=test-$(date '+%d-%m--%H-%M').txt;
BUILDLOG_NAME=log.txt

DIR=/home/sam/datasig/libalgebra_tests
echo "Moving to ${DIR}/build"
cd ${DIR}/build


echo "Starting build"
cmake -DCMAKE_BUILD_TYPE=Release "$@" ..
cmake --build . -- -j6 

echo "Done build"

echo "Running tests"
./test 2>&1 | tee $FNAME;

echo ""
echo "Results saved into $FNAME"

echo ""
echo "Running benchmarks"
./benchmarks



echo "Cleaning up"
rm "test"
rm "benchmarks"



echo "Done"

