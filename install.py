#!/usr/bin/env python3

import os, platform, shutil, subprocess, sys

os.environ["CMAKE_POLICY_VERSION_MINIMUM"] = "3.5"

include_path = ""
lib_path = ""

def is_ubuntu():
    if platform.system() != "Linux":
        return False
    try:
        with open("/etc/os-release") as f:
            return "ubuntu" in f.read().lower()
    except FileNotFoundError:
        return False

def build_optimal_bwt():
    print("###### Configuring and compiling optimalBWT...")
    subprocess.call("make -C external/optimalBWT/".split())

def build_sdsl_lite():
    print("###### Compiling sdsl-lite...")
    shutil.copy("configuration/sdsl-lite/CMakeLists.txt", "external/sdsl-lite/")
    shutil.copy("configuration/sdsl-lite/build.sh", "external/sdsl-lite/build/")
    subprocess.call([
        "external/sdsl-lite/install.sh",
        "external/sdsl-lite/install_dir/",
    ])

def build_remap():
    print("###### Compiling string collection ordering script...")
    include_path = "external/sdsl-lite/install_dir/include"
    lib_path = "external/sdsl-lite/install_dir/lib"
    subprocess.call(
        f"g++ -o remap -std=c++11 -O3 -I{include_path} -L{lib_path} remap.cpp -lsdsl -ldivsufsort -ldivsufsort64".split()
    )

def build_bcr_lcp_gsa():
    print("###### Configuring and Compiling BCR_LCP_GSA...")
    shutil.copy("configuration/Parameters.h", "external/BCR_LCP_GSA/")
    shutil.copy("configuration/makefile", "external/BCR_LCP_GSA/")
    subprocess.call("make DA=1 -C external/BCR_LCP_GSA/".split())

def build_r_index():
    print("###### Configuring and Compiling the r-index...")
    shutil.copy("configuration/r-index/ri-build.cpp", "external/r-index/")
    shutil.copy("configuration/r-index/CMakeLists.txt", "external/r-index/")
    shutil.copy(
        "configuration/r-index/r_index.hpp", "external/r-index/internal/"
    )
    shutil.copy(
        "configuration/r-index/rle_string.hpp", "external/r-index/internal/"
    )

    os.makedirs("external/r-index/build", exist_ok=True)
    subprocess.call(["cmake", ".."], cwd="external/r-index/build")
    subprocess.call(["make"], cwd="external/r-index/build")

def build_big_bwt():
    if is_ubuntu():
        print("###### Compiling Big-BWT (Linux detected)...")
        subprocess.call("make -C external/Big-BWT/".split())
    else:
        print("###### Skipping Big-BWT (Not running on MacOS)...")

# --- Main --- #

def main():
    build_optimal_bwt()
    build_sdsl_lite()
    build_remap()
    build_bcr_lcp_gsa()
    build_r_index()
    build_big_bwt()

##########################
if __name__ == '__main__':
    main()
