#!/PATHCONDAENV/envs/nap/bin/python

import os
import shutil
import sys
import re
#import multiprocessing

if not bool(re.search("chunk", sys.argv[1])):
    exit()

with open(sys.argv[1]) as f:
    lines = f.read().splitlines() 

os.environ['PATH'] += ':/PATHCONDAENV/envs/nap/bin/'
os.environ['JAVA_HOME'] = '/PATHCONDAENV/envs/nap/jre'
os.environ['LD_LIBRARY_PATH'] = '/PATHCONDAENV/envs/nap/jre/lib/amd64/server'

def do_line(line):
    os.system('/PATHTOOLSFOLDER/nap_ccms2/Snap/consensus.R -f ' + 'fusion_res/' + line)

if __name__ == "__main__":
    #pool = multiprocessing.Pool(multiprocessing.cpu_count()-1)
    #pool = multiprocessing.Pool()
    #results = pool.map(do_line, lines)
    for line in lines:
       do_line(line)

