#!/usr/bin/env python3

import sys, subprocess, os

print("###### Compiling remapping script...")
subprocess.call("g++ -o remap -std=c++11 remap.cpp -ldivsufsort -lsdsl".split())

print("###### Compiling BCR_LCP_GSA...")
subprocess.call("cp configuration/Parameters.h external/BCR_LCP_GSA/".split())
subprocess.call("make DA=1 -C external/BCR_LCP_GSA/".split())

print("###### Compiling r-index...")
subprocess.call("cp configuration/ri-build.cpp external/r-index/".split())
subprocess.call("cp configuration/CMakeLists.txt external/r-index/".split())
subprocess.call("cp configuration/r_index.hpp external/r-index/internal/".split())
subprocess.call("cp configuration/rle_string.hpp external/r-index/internal/".split())

subprocess.call("mkdir external/r-index/build".split())
os.chdir("external/r-index/build")
subprocess.call("cmake ..".split())
subprocess.call("make".split())
os.chdir("../../../")

print("###### Compiling optimalBWT...")
subprocess.call("cp configuration/optimalBWT.py external/optimalBWT/".split())
subprocess.call("make -C external/optimalBWT/".split())

print("###### Compiling Big-BWT...")
subprocess.call("make -C external/Big-BWT/".split())

