#!/PATHCONDAENV/envs/nap/bin/python

import os
import shutil
import sys
import re
#import multiprocessing
import math 

with open(sys.argv[1]) as f:
    lines = f.read().splitlines() 

os.environ['PATH'] += ':/PATHCONDAENV/envs/nap/bin/'
os.environ['JAVA_HOME'] = '/PATHCONDAENV/envs/nap/jre'
os.environ['LD_LIBRARY_PATH'] = '/PATHCONDAENV/envs/nap/jre/lib/amd64/server'

def do_line(line):
	os.system("/PATHTOOLSFOLDER/nap_ccms2/Snap/get.metfrag.single.R -f " + 'split_data/' + line + " -O " + "fragmenter_res/" +line)

if __name__ == "__main__":
    #pool = multiprocessing.Pool(math.ceil(multiprocessing.cpu_count()/2))
    #pool = multiprocessing.Pool()
    #results = pool.map(do_line, lines)
    for line in lines:
       do_line(line)

 
