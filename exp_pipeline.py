#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, os, glob, math
from pathlib import Path

Description = """
Script implementing the optimal r-index experimental pipeline
"""

def remove_files(file_list):
    for file_path in file_list:
        Path(file_path).unlink(missing_ok=True)

base_path = os.path.dirname(os.path.abspath(__file__))

BCR_LCP_GSA_exe = base_path + "/external/BCR_LCP_GSA/BCR_LCP_GSA"
r_index_exe = base_path + "/external/r-index/build/ri-build"
r_index_space_exe = base_path + "/external/r-index/build/ri-space"
optbwt_exe = base_path + "/external/optimalBWT/optimalBWT.py"
bigbwt_exe = base_path + "/external/Big-BWT/bigbwt"
remap_exe = base_path + "/remap"

BWT_variants = {"multi": True, "opt": True, "concat": False}

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input_folder', help='folder containing FASTA files', type=str)
    parser.add_argument('output_file',  help='csv output file path', type=str)
    args = parser.parse_args()

    logfile = open(args.output_file+".log","a")

    with open(args.output_file,"w+") as res:

        header = "index_type,dataset_name,dataset_size,no_seq,index_size\n"
        res.write(header)

        dataset_list = os.listdir(args.input_folder)

        for D in dataset_list:

            if D.split(".")[-1] != "fasta":
                continue

            print("############# Processing ", D)

            args.input = args.input_folder + "/" + D
            command = "grep -v >"
            with open(args.input,"rb") as a:
                ps = subprocess.Popen(command.split(), stdin=a, stdout=subprocess.PIPE)
            command = "wc -c"
            output = subprocess.check_output(command.split(), stdin=ps.stdout)
            data_size = str(int(output))

            command = "wc -l"
            with open(args.input,"rb") as a:
                output = subprocess.check_output(command.split(), stdin=a)
            no_seq = int(output)
            if no_seq % 2 != 0:
                no_seq += 1
            no_seq /= 2
            no_seq = int(no_seq)

            if BWT_variants["multi"]:
                #print("Size of the input-order r-index: ",end="")
                command = "{exe} {input} {output}".format(exe = BCR_LCP_GSA_exe,
                                                          input = args.input,
                                                          output= args.input)
                if not execute_command(command,logfile,args.output_file+".log"):
                    return

                command = "{exe} -multidollar {input}".format(exe = r_index_exe,
                                                              input = args.input)
                process = subprocess.check_output(command.split())
                process = str(process).split("\\n")

                res.write("multi-index,"+D+","+data_size+","+str(no_seq)+","+str(process[-2].split(" ")[-1])+"\n")

                remove_files([
                    f"{args.input}.ebwt",
                    f"{args.input}.info",
                    f"{args.input}.len",
                    f"{args.input}.posSA",
                    f"{args.input}.ri",
                ])

            if BWT_variants["opt"]:
                #print("Size of the optimal BWT r-index: ",end="")
                command = "python3 {exe} -a bcr --fasta -b 1000 {input} {output}".format(exe = optbwt_exe, 
                                                                                         input = args.input,
                                                                                         output = args.input)
                if not execute_command(command,logfile,args.output_file+".log"):
                    return

                command = "{exe} {input} {output}".format(exe = remap_exe,
                                                          input = args.input+".optbwt",
                                                          output = args.input+"_remap.fasta")
                if not execute_command(command,logfile,args.output_file+".log"):
                    return

                command = "{exe} {input} {output}".format(exe = BCR_LCP_GSA_exe,
                                                          input = args.input+"_remap.fasta",
                                                          output= args.input+"_remap.fasta")
                if not execute_command(command,logfile,args.output_file+".log"):
                    return

                command = "{exe} -multidollar {input}".format(exe = r_index_exe,
                                                              input = args.input+"_remap.fasta")
                process = subprocess.check_output(command.split())
                process = str(process).split("\\n")

                res.write("opt-index,"+D+","+data_size+","+str(no_seq)+","+str(process[-2].split(" ")[-1])+"\n")

                remove_files([
                    f"{args.input}.log",
                    f"{args.input}.info",
                    f"{args.input}.len",
                    f"{args.input}.optbwt",
                    f"{args.input}._remap.fasta",
                    f"{args.input}._remap.fasta.ebwt",
                    f"{args.input}._remap.fasta.info",
                    f"{args.input}._remap.fasta.len",
                    f"{args.input}._remap.fasta.posSA",
                    f"{args.input}._remap.fasta.ri",
                ])

            if BWT_variants["concat"]:
                #print("Size of the original r-index: ",end="")
                command = "grep -v >"
                with open(args.input, 'rb', 0) as a, open(args.input+".flat", 'w') as b:
                    rc = subprocess.call(command.split(), stdin=a, stdout=b)

                command = "python3 {exe} -w 10 -p 50 -s -e {input}".format(exe = bigbwt_exe,
                                                                           input = args.input+".flat")
                if not execute_command(command,logfile,args.output_file+".log"):
                    return

                command = "{exe} {input}".format(exe = r_index_exe,
                                                 input = args.input+".flat")
                process = subprocess.check_output(command.split())
                process = str(process).split("\\n")

                res.write("concat-index,"+D+","+data_size+","+str(no_seq)+","+str(process[-2].split(" ")[-1])+"\n")

                remove_files([
                    f"{args.input}.flat",
                    f"{args.input}.flat.ri",
                    f"{args.input}.flat.ssa",
                    f"{args.input}.flat.log",
                    f"{args.input}.flat.esa",
                    f"{args.input}.flat.bwt",
                ])

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