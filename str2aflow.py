#!/usr/bin/python

f=open('struct_enum.out','r')
g=open('struct_enum.aflow','w')
# Skip the first 12 lines of the struct_enum.out file
skip = range(12)
for i in skip:
    s = f.readline()
    g.write(s.replace('\n',' \\n\\\n'))
for s in f:
    g.write(s.replace('\n',' ??? \\n\\\n'))
f.close
g.close


