#!/usr/bin/env python3

import os, platform, shutil, subprocess, sys

os.environ["CMAKE_POLICY_VERSION_MINIMUM"] = "3.5"

def is_linux():
    return platform.system() == "Linux"

def run(command, cwd=None):
    if subprocess.call(command, cwd=cwd) != 0:
        sys.exit("###### Error: command failed: " + " ".join(command))

def build_optimal_bwt():
    print("###### Configuring and compiling optimalBWT...")
    shutil.copy("configuration/optimalBWT/optimalBWT.py", "external/optimalBWT/")
    shutil.copy("configuration/BCR_LCP_GSA/makefile", "external/optimalBWT/external/BCR_LCP_GSA/")
    run("make -C external/optimalBWT/".split())

def build_sdsl_lite():
    print("###### Compiling sdsl-lite...")
    shutil.copy("configuration/sdsl-lite/CMakeLists.txt", "external/sdsl-lite/")
    run([
        "external/sdsl-lite/install.sh",
        "external/sdsl-lite/install_dir/",
    ])

def build_remap():
    print("###### Compiling string collection ordering script...")
    include_path = "external/sdsl-lite/install_dir/include"
    lib_path = "external/sdsl-lite/install_dir/lib"
    run(
        f"g++ -o remap -std=c++11 -O3 -I{include_path} -L{lib_path} remap.cpp -lsdsl -ldivsufsort -ldivsufsort64".split()
    )

def build_bcr_lcp_gsa():
    print("###### Configuring and Compiling BCR_LCP_GSA...")
    shutil.copy("configuration/Parameters.h", "external/BCR_LCP_GSA/")
    shutil.copy("configuration/BCR_LCP_GSA/makefile", "external/BCR_LCP_GSA/")
    run("make DA=1 -C external/BCR_LCP_GSA/".split())

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
    run(["cmake", ".."], cwd="external/r-index/build")
    run(["make"], cwd="external/r-index/build")

def build_big_bwt():
    if is_linux():
        print("###### Compiling Big-BWT (Linux detected)...")
        run("make -C external/Big-BWT/".split())
    else:
        print("###### Skipping Big-BWT (only supported on Linux; needed only by exp_pipeline.py)...")

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
