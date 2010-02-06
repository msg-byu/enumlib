#!/usr/bin/python

import os,unclefuncs,sys
from math import sqrt

#print 'You will need to make sure the file makestr.2d is in this folder.  It is obtained by compiling makestr2d.f90'
#verify = str(raw_input('Is it there?'))
#if verify != 'yes':
#    sys.exit('Do it then!')
#if os.path.isfile('inputvasp.txt') == False:
#    print 'You are lacking an input file called inputvasp.txt, which contains a list of the structures for which you desire to compute energies'
#    sys.exit()
input = unclefuncs.readfile('inputvasp.txt')

poscars = input[0].split()

if '-' in input[0].split():
    poscars = []
    for i in range(int(input[0].split()[0]),int(input[0].split()[2])+1):
        poscars.append(str(i))


struct_enum = unclefuncs.readfile('struct_enum.out')

first = int(struct_enum[15].split()[0])
#print poscars
#if first != 1:
#    for i in range(len(poscars)):
#        poscars[i] = str(int(poscars[i]) + (first - 1))
#print poscars

writelines = []
for j in poscars:
    os.system('./makestr.2d struct_enum.out ' + j)
    
    structnum = str(int(j) + first - 1)
    print structnum
    if int(structnum) > 9999:
        pos = unclefuncs.readfile('vasp.0'+ structnum)
    elif int(structnum) > 999:
        pos = unclefuncs.readfile('vasp.00'+ structnum)
    elif int(structnum) > 99:
        pos = unclefuncs.readfile('vasp.000'+ structnum)
    elif int(structnum) > 9:
        pos = unclefuncs.readfile('vasp.0000'+ structnum)
    else:
        pos = unclefuncs.readfile('vasp.00000'+ structnum)

    latvecone = pos[3].split()
    latvectwo = pos[4].split()
    upperlimitx = float(latvecone[1]) + float(latvectwo[1])
    upperlimity = float(latvecone[2]) + float(latvectwo[2])
    if upperlimitx > upperlimity:
        upperlimit = int(upperlimitx)
    else: 
        upperlimit = int(upperlimity)
    latvecs = [latvecone[1] + ' ' + latvecone[2],latvectwo[1] + ' ' + latvectwo[2]]

    numofbasis = pos[5].split()
    if len(numofbasis) == 1:
        numofatoms = int(numofbasis[0])
    elif len(numofbasis) == 2:
        numofatoms = int(numofbasis[0]) + int(numofbasis[1])

    aatoms = []
    batoms = []
    atoms = []
    for k in range(7,7+numofatoms):
        atoms.append([float(pos[k].split()[1]),float(pos[k].split()[2])])

    for k in range(7,7+int(numofbasis[0])):
        aatoms.append([float(pos[k].split()[1]),float(pos[k].split()[2])])
    for k in range(7+int(numofbasis[0]),7+numofatoms):
        batoms.append([float(pos[k].split()[1]),float(pos[k].split()[2])])


    allatoms = []

    for l in atoms:
        for i in range(0,int(upperlimit) + 1):
            for z in range(0,int(upperlimit) + 1):
                shiftvec = [i, z]
                newatom = [float(l[0]) + shiftvec[0], float(l[1]) + shiftvec[1]]
                allatoms.append(newatom)


    finalaatomslist = []
    finalbatomslist = []
    for m in allatoms:
        mapped = unclefuncs.mapback2d(latvecs,m,atoms,upperlimit)

        if mapped in aatoms:
            finalaatomslist.append(m)

        elif mapped in batoms:
            finalbatomslist.append(m)

    xlista = []
    xlistb = []

    ylista = []
    ylistb = []

    for f in finalaatomslist:
        xlista.append(f[0])
        ylista.append(f[1])
    for g in finalbatomslist:
        xlistb.append(g[0])
        ylistb.append(g[1])
        

    status = 'continue'
#    for h in xlista:
#        if h not in xlistb:
#            status = 'done'
#            break

#    if status == 'continue':
#        for r in ylista:
#            if r not in ylistb:
#                status = 'done'

    if status == 'continue':
        for p in xlistb:
            if p not in xlista:
                status = 'done'

    if status == 'continue':
        for q in ylistb:
            if q not in ylista:
               status = 'done'
    if status == 'done':
        energy = '  1\n'
    else:
        energy = '  0\n'

    writelines.append(structnum + energy)

os.system('rm vasp.*')

unclefuncs.writefile('gssreal.out',writelines)
