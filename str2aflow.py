#!/usr/bin/python
# This script takes a struct_enum.out file and makes the format easy 
# for aflow to slurp up and use
import os, sys
import re

tag = 'enum'
r1 = re.compile('^(\s+)([0-9]+.*)')
f=open('struct_enum.out','r')
#f=open('strtest')
print 'Making "struct_enum.aflow" file...'

# Read the system type from the first line
systype = f.readline()
systype = str(systype[0:3]) # Keep just the first three characters
systype = re.sub(' ','',systype) # Strip trailing blanks
g=open('struct_enum.aflow.'+systype,'w')

# Skip the next 11 lines of the struct_enum.out file
skip = range(11)
for i in skip:
    s = f.readline()
    g.write(s)
for s in f:
    m = r1.match(s)
    g.write(re.sub(r1,'['+tag+systype+']'+m.group(1)+tag+systype+m.group(2)+' ???',s))
f.close
g.close
#print 'Zipping up struct_enum.aflow file...'
#os.system('gzip -9 struct_enum.aflow.'+systype)
