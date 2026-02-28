#!/bin/bash
mkdir -p build && cd build
cmake ..
make -j4
./test_harness
./test_cdt
