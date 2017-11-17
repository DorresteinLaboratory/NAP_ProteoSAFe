#!/PATHCONDAENV/envs/nap/bin/python

import os
import shutil
import sys
import re
import math 

def main():
    os.system("/PATHTOOLSFOLDER/nap_ccms2/Snap/load.frag.R -d " + sys.argv[1] + " -f " + sys.argv[2] + " -O " + sys.argv[3])

    fls = os.listdir('fragmenter_res')
    fls = [x for x in fls  if bool(re.match('line', x))]
    fls = sorted(fls, key=lambda e: int(re.sub('\\D', '', e)))

    chunksize = 15
    if len(fls)/15 > 30: 
        chunksize = math.ceil(len(fls)/30) 
    
    k = 1
    for i in range(0, len(fls), chunksize):
        to_save = fls[i:i+chunksize]
        with open('chunk_res/chunk' + str(k) + '.txt', 'w') as fp:
            for item in to_save:
                fp.write("%s\n" % item)
     
        k += 1

if __name__ == "__main__":
    main()

