#!/PATHCONDAENV/envs/nap/bin/python

import os
import shutil
import sys


def main():
    input_par = sys.argv[1]

    os.environ["JAVA_HOME"] = '/PATHCONDAENV/envs/nap/jre'
    os.environ["LD_LIBRARY_PATH"] = '/PATHCONDAENV/envs/nap/jre/lib/amd64/server'
    if sys.argv[2] != "" and sys.argv[3] != "":
        input_filename = sys.argv[2]
        input_filename2 = sys.argv[3]
        os.system('/PATHTOOLSFOLDER/nap_ccms2/split_network.R '  + input_par + ' ' + input_filename + ' ' + input_filename2)
    elif sys.argv[2] != "":
        input_filename = sys.argv[2]
        os.system('/PATHTOOLSFOLDER/nap_ccms2/split_network.R '  + input_par + ' ' + input_filename)
    else:
        os.system('/PATHTOOLSFOLDER/nap_ccms2/split_network.R ' + input_par + ' ' + 'split_data')

if __name__ == "__main__":
    main()

