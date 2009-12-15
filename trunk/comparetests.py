#!/usr/bin/python

import os, subprocess
import glob
from os import system
from shutil import copyfile
import re
import string
import sys

system('clear') # Clear the screen and start fresh
rs=system('make')  # need the exit status here...
if rs!=0: sys.exit("\n\n *** Compilation of enumlib failed *** \n")
rs=system('make makestr.x')  # need the exit status here...
if rs!=0: sys.exit("\n\n *** Compilation of makestr.x failed *** \n")

#tests=glob.glob('tests/struct_enum.in.*')
#f=open('tests/list')

system('rm -r vasp.0007')
print "\n *** testing on of Rod's gang of 13 ***"
system('./makestr.x ./tests/13struct.out 7')
rs=system('diff -q vasp.0007 tests/vasp.0007_from2x3_Rods_case')
if rs!=0:
    sys.exit("\n --- Failure in the first of Rod's 2x3 cases ---\n")
else:
    print "\n <<< First test of makestr.x passed >>>\n"

# Loop over the 17 fcc superstructures that have 2-4 atoms/cell. For pictures, see
# "Where are nature's missing structures" Nat. Mat. 2008 or "Algorithm for generating
# derivative superstructures" PRB 2008.
poscar = glob.glob('tests/*.poscar')
f = open('tests/17list')
rs=system('make find_structure_in_list.x')  # need the exit status here...
print "\n *** Testing compare code on first 17 fcc structures ***"
if rs!=0: sys.exit("\n\n *** Compilation of find_structure_in_list.x failed *** \n")
for ip in poscar:
    #print ip
    command = './find_structure_in_list.x '+ip+' tests/struct17fcc.out'
    #process = subprocess.Popen('pwd')
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    os.waitpid(process.pid, 0)
    #print process.stdout.read()
    #print process.stdout.read().strip() \w+([0-9]+)
    r1 = re.compile('Structure #: \s+([0-9]+)')
    m = r1.match(process.stdout.read().strip())
    nStr =  m.group(1)
    #print f.readline().strip()+"X"+nStr+"X
    if(f.readline().strip()!=nStr):
        sys.exit("\n\nFailure in test")
    #print "test passed"

print "\n <<< Test of 17 fcc with compare code passed >>>\n"
