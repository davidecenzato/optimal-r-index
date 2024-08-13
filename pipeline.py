#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, os, glob, math

Description = """
Tool for the optimalBWT journal
"""

base_path = os.path.dirname(os.path.abspath(__file__))

BCR_LCP_GSA_exe = base_path + "/external/BCR_LCP_GSA/BCR_LCP_GSA"
r_index_exe = base_path + "/external/r-index/build/ri-build"
r_index_space_exe = base_path + "/external/r-index/build/ri-space"
optbwt_exe = base_path + "/external/optimalBWT/optimalBWT.py"
bigbwt_exe = base_path + "/external/Big-BWT/bigbwt"
remap_exe = base_path + "/remap"

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input_folder', nargs='?', help='a folder containing fasta files', type=str)
    parser.add_argument('output_file', nargs='?', help='csv output file name', type=str)
    parser.add_argument('--multi',  help='compute multidollarBWT r-index size (def. False)',action='store_true')
    parser.add_argument('--opt',  help='compute optimalBWT r-index size (def. False)',action='store_true')
    parser.add_argument('--concat',  help='compute concatBWT r-index size (def. False)',action='store_true')
    args = parser.parse_args()

    logfile = open(args.output_file+".log","a")

    with open(args.output_file,"w+") as res:

        header = "index_type,dataset_name,dataset_size,no_seq,index_size\n"
        res.write(header)

        dataset_list = os.listdir(args.input_folder)

        for D in dataset_list:

            if D.split(".")[-1] != "fasta":
                continue

            print("############# Processing ",D)

            args.input = args.input_folder + "/" + D
            command = "grep -v >"
            with open(args.input,"rb") as a:
                ps = subprocess.Popen(command.split(), stdin=a, stdout=subprocess.PIPE)
            command = "wc -c"
            output = subprocess.check_output(command.split(), stdin=ps.stdout)
            data_size = str(int(output))
            #print("DATA SIZE=",data_size)

            command = "wc -l"
            with open(args.input,"rb") as a:
                output = subprocess.check_output(command.split(), stdin=a)
            no_seq = int(output)
            if no_seq % 2 != 0:
                no_seq += 1
            no_seq /= 2
            no_seq = int(no_seq)
            #print("NO_SEQ=",no_seq)

            if args.multi:
                print("Size of the multidollar BWT r-index: ",end="")
                command = "{exe} {input} {output}".format(exe = BCR_LCP_GSA_exe, input = args.input, output= args.input)
                #print("###########",command)
                if not execute_command(command,logfile,args.output_file+".log"):
                    return

                command = "{exe} -multidollar {input}".format(exe = r_index_exe, input = args.input)
                #print("###########",command)
                process = subprocess.check_output(command.split())
                process = str(process)
                process = process.split("\\n")
                #print("Size=",process[-2].split(" ")[-1])
                print(process[-2].split(" ")[-1])

                res.write("multi-index,"+D+","+data_size+","+str(no_seq)+","+str(process[-2].split(" ")[-1])+"\n")

                os.remove(args.input+".ebwt")
                os.remove(args.input+".info")
                os.remove(args.input+".len")
                os.remove(args.input+".posSA")
                os.remove(args.input+".ri")
                os.remove(args.input+".da")

            if args.opt:
                print("Size of the optimal BWT r-index: ",end="")
                command = "python3 {exe} -a bcr --fasta -b 1000 {input} {output}".format(exe = optbwt_exe, input = args.input, output = args.input)
                #print(command)
                if not execute_command(command,logfile,args.output_file+".log"):
                    return

                command = "{exe} {input} {output}".format(exe = remap_exe, input = args.input+".optbwt", output = args.input+"_remap.fasta")
                #print(command)
                if not execute_command(command,logfile,args.output_file+".log"):
                    return

                command = "{exe} {input} {output}".format(exe = BCR_LCP_GSA_exe, input = args.input+"_remap.fasta", output= args.input+"_remap.fasta")
                #print(command)
                if not execute_command(command,logfile,args.output_file+".log"):
                    return

                command = "{exe} -multidollar {input}".format(exe = r_index_exe, input = args.input+"_remap.fasta")
                #print("###########",command)
                process = subprocess.check_output(command.split())
                process = str(process)
                process = process.split("\\n")
                #print("Size=",process[-2].split(" ")[-1])
                print(process[-2].split(" ")[-1])

                res.write("opt-index,"+D+","+data_size+","+str(no_seq)+","+str(process[-2].split(" ")[-1])+"\n")

                os.remove(args.input+".log")
                os.remove(args.input+".info")
                os.remove(args.input+".len")
                os.remove(args.input+".optbwt")
                os.remove(args.input+"_remap.fasta")
                os.remove(args.input+"_remap.fasta.ebwt")
                os.remove(args.input+"_remap.fasta.info")
                os.remove(args.input+"_remap.fasta.len")
                os.remove(args.input+"_remap.fasta.posSA")
                os.remove(args.input+"_remap.fasta.ri")
                os.remove(args.input+"_remap.fasta.da")

            if args.concat:
                print("Size of the original r-index: ",end="")
                command = "grep -v >"
                with open(args.input, 'rb', 0) as a, open(args.input+".flat", 'w') as b:
                    rc = subprocess.call(command.split(), stdin=a, stdout=b)

                command = "python3 {exe} -w 10 -p 50 -s -e {input}".format(exe = bigbwt_exe, input = args.input+".flat")
                if not execute_command(command,logfile,args.output_file+".log"):
                    return

                command = "{exe} {input}".format(exe = r_index_exe, input = args.input+".flat")
                #print("###########",command)
                process = subprocess.check_output(command.split())
                process = str(process)
                process = process.split("\\n")
                #print("Size=",process[-2].split(" ")[-1])
                print(process[-2].split(" ")[-1])

                res.write("concat-index,"+D+","+data_size+","+str(no_seq)+","+str(process[-2].split(" ")[-1])+"\n")

                os.remove(args.input+".flat")
                os.remove(args.input+".flat.ri")
                os.remove(args.input+".flat.ssa")
                os.remove(args.input+".flat.log")
                os.remove(args.input+".flat.esa")
                os.remove(args.input+".flat.bwt")

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

if __name__ == '__main__':
    main()