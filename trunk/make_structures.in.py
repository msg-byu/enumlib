#!/usr/bin/python

import os, subprocess
import glob
from os import system
import string
import sys

# Read in a struct_enum-format file with structures defined. These will be converted to the
# vasp-like format that UNCLE uses
if len( sys.argv) < 2:
    sys.exit("Script needs an argument---file containing structures")
#strfile = open(sys.argv[1])
for n in range(1,1300):
    rs=system('./makestr.x '+sys.argv[1]+' '+str(n)+' >& /dev/null')
    if rs!=0:
         break
print "\n Read in a total of "+str(n-1)+" stuctures from "+sys.argv[1]
system("perl -pi -e 's/scale factor/1.0/g' vasp.*")
unclefile = open("structures.in",'w')
unclefile.write("scalar\nperatom\nnoweights\n")
vnum =glob.glob('vasp.*')
for i in vnum:
    f = open(i)
    str = f.read()
    unclefile.write(str)
    unclefile.write("0 0 \n #\n")
unclefile.close
system("rm vasp.*")
print "\n Structures concatenated into structures.in\n"
