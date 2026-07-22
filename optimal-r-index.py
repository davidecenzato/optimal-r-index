#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, os, glob, math

Description = """
Tool for computing the optimal r-index of a string collection
"""

base_path = os.path.dirname(os.path.abspath(__file__))

BCR_LCP_GSA_exe = base_path + "/external/BCR_LCP_GSA/BCR_LCP_GSA"
r_index_exe = base_path + "/external/r-index/build/ri-build"
optbwt_exe = base_path + "/external/optimalBWT/optimalBWT.py"
remap_exe = base_path + "/remap"

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input',  nargs='?', help='input FASTA path', type=str)
    parser.add_argument('output', nargs='?', help='output file path', type=str)
    parser.add_argument('--algo', help='algorithm for computing the optimal BWT (sais|bcr). Default: sais',
                                  default="sais", type=str, choices=["sais", "bcr"])
    parser.add_argument('--keep', help='keep temporary files',action='store_true')
    args = parser.parse_args()

    logfile = open(args.input+".log","a")

    command = "python3 {exe} -a {algo} --fasta -b 1000 {input} {output}".format(exe = optbwt_exe,
                                                                                algo = args.algo,
                                                                                input = args.input,
                                                                                output = args.input)
    if not execute_command(command,logfile,args.input+".log"):
        return

    command = "{exe} {input} {output}".format(exe = remap_exe,
                                              input = args.input+".optbwt",
                                              output = args.input+"_sorted.fasta")
    if not execute_command(command,logfile,args.input+".log"):
        return

    command = "{exe} {input} {output}".format(exe = BCR_LCP_GSA_exe,
                                              input = args.input+"_sorted.fasta",
                                              output= args.output)
    if not execute_command(command,logfile,args.input+".log"):
        return

    command = "{exe} -multidollar {input}".format(exe = r_index_exe,
                                                  input = args.output)
    if not execute_command(command,logfile,args.input+".log"):
        return

    if not args.keep:
        files_to_remove = [
            args.input + ".optbwt",
            args.input + "_sorted.fasta",
            *(f"{args.output}.{ext}" for ext in ["ebwt", "info", "len", "posSA"]),
        ]

        for filepath in files_to_remove:
            os.remove(filepath)

    logfile.close()

# execute command: return True is everything OK, False otherwise
def execute_command(command,logfile=None,logfile_name=None,env=None):
    try:
        #subprocess.run(command.split(),stdout=logfile,stderr=logfile,check=True,env=env)
        subprocess.check_call(command.split(),stdout=logfile,stderr=logfile,env=env)
    except subprocess.CalledProcessError:
        print("Error executing command line:")
        print("\t"+ command)
        print("Check log file: " + logfile_name)
        return False

    return True

##########################
if __name__ == '__main__':
    main()