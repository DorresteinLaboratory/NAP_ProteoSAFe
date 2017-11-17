#!/PATHCONDAENV/envs/nap/bin/python

import os
import shutil
import json
import sys

def main():
    os.environ['PATH'] += ':/PATHCONDAENV/envs/nap/bin/'
    os.environ['JAVA_HOME'] = '/PATHCONDAENV/envs/nap/jre'
    os.environ['LD_LIBRARY_PATH'] = '/PATHCONDAENV/envs/nap/jre/lib/amd64/server'
    os.system('/PATHTOOLSFOLDER/nap_ccms2/Snap/metfrag.out.R')

if __name__ == "__main__":
    main()

