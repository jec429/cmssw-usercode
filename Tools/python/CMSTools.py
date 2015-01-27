#!/usr/bin/env python

import os,glob,subprocess
def file_list(name,local):
    if local:
        fs = glob.glob(name)
        myfilelist = []
        for f in fs:
        # add as many files as you wish this way
            myfilelist.append(f.replace('/eos/uscms',''))
        return myfilelist
    else:
        #print name
        myfilelist = []
        os.system('das_client.py --query="file dataset=' + name + '" --limit=0 > pappy.txt')
        f = open('pappy.txt')
        for line in f:
            myfilelist.append(line)        
        os.system('rm pappy.txt')
        return myfilelist
