#!/usr/bin/python

import os,unclefuncs,sys
from math import sqrt

print 'You will need to copy the file makestr.2d to your bin for this program to work'
verify = str(raw_input('Is it there?'))
if verify != 'yes':
    sys.exit('Do it then!')
if os.path.isfile('inputvasp.txt') == False:
    print 'You are lacking an input file called inputvasp.txt, which contains a list of the structures for which you desire to compute energies'
    sys.exit()
input = unclefuncs.readfile('inputvasp.txt')

poscars = input[0].split()

if '-' in input[0].split():
    poscars = []
    for i in range(int(input[0].split()[0]),int(input[0].split()[2])+1):
        poscars.append(str(i))

writelines = []
for j in poscars:
    os.system('makestr.test struct_enum.out ' + j)
    
    if int(j) > 9999:
        pos = unclefuncs.readfile('vasp.0'+ j)
    elif int(j) > 999:
        pos = unclefuncs.readfile('vasp.00'+ j)
    elif int(j) > 99:
        pos = unclefuncs.readfile('vasp.000'+ j)
    elif int(j) > 9:
        pos = unclefuncs.readfile('vasp.0000'+ j)
    else:
        pos = unclefuncs.readfile('vasp.00000'+ j)

    latvecone = pos[3].split()
    latvectwo = pos[4].split()
    upperlimitx = float(latvecone[1]) + float(latvectwo[1])
    upperlimity = float(latvecone[2]) + float(latvectwo[2])
    if upperlimitx > upperlimity:
        upperlimit = upperlimitx
    else: 
        upperlimit = upperlimity
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
        for i in range(0,upperlimit + 1):
            for z in range(0,upperlimit + 1):
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
    for h in xlista:
        if h not in xlistb:
            status = 'done'
            break

    if status == 'continue':
        for r in ylista:
            if r not in ylistb:
                status = 'done'

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

    writelines.append(j + energy)

os.system('rm vasp.*')

unclefuncs.writefile('gssreal.out',writelines)
