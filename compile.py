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
subprocess.call("cp configuration/r-index.hpp external/r-index/internal/".split())
subprocess.call("cp configuration/rle_string.hpp external/r-index/internal/".split())

subprocess.call("mkdir external/r-index/external".split())
os.chdir("mkdir external/r-index/external")
subprocess.call("git clone https://github.com/simongog/sdsl-lite".split())
os.chdir("sdsl-lite")
subprocess.call("./install installed".split())
os.chdir("../../../../")
subprocess.call("mkdir external/r-index/build".split())
os.chdir("external/r-index/build")
subprocess.call("cmake ..".split())
subprocess.call("make".split())
os.chdir("../../../")

print("###### Compiling optimalBWT...")
subprocess.call("make -C external/optimalBWT/".split())

