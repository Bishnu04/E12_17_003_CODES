#!/apps/python/3.4.3/bin/python3

# Toshiyuki Gogami
# Nov 2, 2018

import sys
import time, os.path
from subprocess import call
#import concurrent.futures
from logging import StreamHandler, Formatter, INFO, getLogger
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures.process import ProcessPoolExecutor 
import numpy as np

trigflag = 5 # (coin)
#trigflag = 4 # (RHRS)
#trigflag = 1 # (LHRS)
nworkers=4 # no of computer for work


# the following is the file that contains input root file information
#runfile = "H_149_220.dat"
#runfile = "H_480_542.dat"
#runfile = "HT_552_716.dat" 
#runfile = "T_221_368.dat"
#runfile = "T_369_830.dat"
runfile = "He_all.dat"
#runfile = "test_list.dat" # for test purpose

thisfile = "nude3.py"

def nude_start(command):
    time.sleep(1.0)
    call(command,shell=True)

def main():
    comlist = []
    #inputfile = open("h2_2.dat","r")
    inputfile = open(runfile,"r")
    lines = inputfile.readlines()
    for line in lines:
        data = line.split()
        #com = "./nude " + data[0]+ " " +data[1]
        com = "root -l -q \"nude3.cc(" + data[0]+ "," +data[1] + "," + str(trigflag)
        com2 = com + ")\";"
        #call(com,shell=True)
        comlist.append(com2)
    with ProcessPoolExecutor(max_workers=nworkers) as executor:
        executor.map(nude_start,comlist)


stime = time.time()
main()
print("\n Jobs were done in %.0f sec \n" % float(time.time()-stime))
