#!/usr/bin/python

import os,random
from math import pi,sqrt
def crossprod(vec1,vec2):
    vector = [vec1[1]*vec2[2]-vec1[2]*vec2[1],-(vec1[0]*vec2[2]-vec1[2]*vec2[0]),vec1[0]*vec2[1]-vec1[1]*vec2[0]]
    return vector

def dotprod(vec1,vec2):
    vector = vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]
    return vector

def magnitude(vector):
    mag=sqrt(vector[0]**2+vector[1]**2+vector[2]**2)
    return mag

def magnitude2d(vector):
    mag=sqrt(vector[0]**2+vector[1]**2)
    return mag


def recip(latvecs):
    crossprod1=crossprod(latvecs[1],latvecs[2])
    crossprod2=crossprod(latvecs[2],latvecs[0])
    crossprod3=crossprod(latvecs[0],latvecs[1])
    dotprod1=dotprod(latvecs[0],crossprod1)
    dotprod2=dotprod(latvecs[1],crossprod2)
    dotprod3=dotprod(latvecs[2],crossprod3)
    component1=[x/dotprod1 for x in crossprod1]
    component2=[x/dotprod2 for x in crossprod2]
    component3=[x/dotprod3 for x in crossprod3]
    recipcell =  [component1,component2,component3]
    return recipcell



def kmesh(recipunitcell,resolution):
    mesh1=[[recipunitcell[x][y]/resolution for y in range(3)] for x in range(3)]
    mesh=[str(mesh1[x][0])+ ' ' + str(mesh1[x][1]) + ' ' + str(mesh1[x][2]) for x in range(3)]
    return mesh


def findenergy():
    if os.path.isfile('OUTCAR') == True:
        outcarlines = readfile('OUTCAR')
    else:
        outcarlines = []

    count = 0
    energylines = [];

    for i in outcarlines:
        list = i.split()
        if 'TOTEN' in list:
            energylines.append(count)
        count = count + 1
    
    if len(energylines) > 0:
        last = energylines[-1]
        energyline = outcarlines[last]
        energy = float(energyline.split()[4])
    else:
        energy = '????????'
    return energy

def formationenthalpy(mixture,a,b,aatoms,batoms,totalatoms):
    totatoms = aatoms + batoms
    e = mixture/totatoms - (a/totalatoms * aatoms/totatoms + b/totalatoms * batoms/totatoms)
    return e
    
def updatestructure(structurenumber,formationenthalpy):
    structuresread = open('structures.in','r')
    number = 0
    count = 0
    for i in structuresread:
        if i == ('formation enthalpy' + str(structurenumber) + '\n'):
            number = count
        count = count + 1
    structuresread.close()
    
    
    readlines = readfile('structures.in')
    if number != 0:
        readlines[number] = str(formationenthalpy) + '\n'

    writefile('structures.in',readlines)
    
def vegardslaw(constone,consttwo,largelatpar,smalllatpar):
    latpar =  largelatpar - (largelatpar - smalllatpar) * float(constone)/(float(constone) + float(consttwo))
    return latpar

def getscore():
    if os.path.isfile('finalcvs.out')==True:
        output = open('finalcvs.out','r')
    else:
        output = []

    count = 0
    found = 0
    for i in output:
        if i == " Final Cross Validation Score:\n":
            found = count + 1
        count = count + 1
    if os.path.isfile('finalcvs.out') == True:
        output.close()
        lines = readfile('finalcvs.out')
    
    if found != 0:
        score = float(lines[found])
    else:
        score = '???????????'
    return score

def conv(file,keyword):
    test = open(file,'r')

    count = 0
    energylines = [];

    for i in test:
        list = i.split()
        if keyword in list:
            energylines.append(count)
        count = count + 1
    test.close()

    test1 = open(file,'r')
    lines = test1.readlines()
    test.close()

    energies = []

    for i in energylines:
        energies.append(lines[i])


#    energy = float(energyline.split()[4])
    return energies


def randomize(numofstructs,path):
    print numofstructs
    
    data = range(2,numofstructs + 2)
    ranlist = []
    print data

    for i in range(1,numofstructs + 1):
        print i
        ran = random.choice(data)
        ranlist.append(ran)
        data.remove(ran)
    
    structlines = readfile(path + '/structures.orig')

    list = structuresindices(structlines)

    structlist = ['peratom\n']

    constonelines = ''.join(structlines[1:list[1]])
    structlist.append(constonelines)
    
    for i in ranlist:
        lines = ''.join(structlines[list[i-1]:list[i]])
        structlist.append(lines)
        print structlist
        print '\n'

    
    consttwolines = ''.join(structlines[list[-1]:])
    structlist.append(consttwolines)
    

    writefile(path +'/structures.in',structlist)


def readfile(file):
    openfile = open(file,'r')
    lines = openfile.readlines()
    openfile.close()
    return lines

def writefile(file,lines):
    openfile = open(file,'w')
    openfile.writelines(lines)
    openfile.close()

def structuresindices(structlines):
    index = 0
    list = []
    for i in structlines:
        if '#-----' in i:
            list.append(index)
        index = index + 1
    return list

def getindices(file,whattolookfor):
    index = 0
    list = []
    lines = readfile(file)
    for i in lines:
        if whattolookfor in i:
            list.append(index)
        index = index + 1
    return list


def placeunclefiles(homedir):
    if os.path.isfile(homedir + '/GApar.in') == True:
        os.system('cp ' + homedir + '/GApar.in .')
    if os.path.isfile(homedir + '/lat.in') == True:
        os.system('cp ' + homedir + '/lat.in .')
    if os.path.isfile(homedir + '/fitpar.in') == True:
        os.system('cp ' + homedir + '/fitpar.in .')
    if os.path.isfile(homedir + '/control.in') == True:
        os.system('cp ' + homedir + '/control.in .')
    if os.path.isfile(homedir + '/bulkpar.in') == True:
        os.system('cp ' + homedir + '/bulkpar.in .')
    if os.path.isfile(homedir + '/groundstatesearch.in') == True:
        os.system('cp ' + homedir + '/groundstatesearch.in .')
    if os.path.isfile(homedir + '/MCpar.in') == True:
        os.system('cp ' + homedir + '/MCpar.in .')
    if os.path.isfile(homedir + '/finalecis.out') == True:
        os.system('cp ' + homedir + '/finalecis.out .')
    if os.path.isfile(homedir + '/figures.out') == True:
        os.system('cp ' + homedir + '/figures.out .')
    if os.path.isfile(homedir + '/finalcvs.out') == True:
        os.system('cp ' + homedir + '/finalcvs.out .')
    if os.path.isfile(homedir + '/struct_enum.out') == True:
        os.system('cp ' + homedir + '/struct_enum.out .')
    if os.path.isfile(homedir + '/genalgsummary.dat') == True:
        os.system('cp ' + homedir + '/genalgsummary.dat .')
    if os.path.isfile(homedir + '/listchildren.dat') == True:
        os.system('cp ' + homedir + '/listchildren.dat .')
    if os.path.isfile(homedir + '/listgeneration.dat') == True:
        os.system('cp ' + homedir + '/listgeneration.dat .')
    if os.path.isfile(homedir + '/population.out') == True:
        os.system('cp ' + homedir + '/population.out .')
    if os.path.isfile(homedir + '/gss.out') == True:
        os.system('cp ' + homedir + '/gss.out .')
    
def concentrations(list):
    uniqueconclist = uniquelist(list,1)
        
    if os.path.isdir('concentrations') == False:
        os.mkdir('concentrations')
    else:
        os.system('rm -rf concentrations')
        os.mkdir('concentrations')
    os.chdir('concentrations')
    
    for i in uniqueconclist:
        specificconcs = []
        for j in list:
            if j.split()[1] == i and len(j.split()) > 2:
                specificconcs.append(j)
        specificconcs = sorted(specificconcs)
        writefile('concentration' + str(i),specificconcs)
    os.chdir('../')

def numberofstructures(list):
    uniqueconclist = uniquelist(list,1)
    os.chdir('concentrations')
    

    for l in uniqueconclist:
        total = 0
        struct = ['structure #      # of occcurences\n--------------------------------------------\n']
        plot = []
        lines = readfile('concentration' + l)
        uniquestructlist = uniquelist(lines,0)
        for j in lines:
            struct.append(j.split()[0])
        
        for i in uniquestructlist:
            number = struct.count(i)
            total = total + number
            plot.append(i.rjust(5) + '   ' + str(number).rjust(3) + '\n')
        plot.append('total --> ' + str(total) + '\n')
        writefile('plot' + l,plot)

def uniquelist(list,index):
    uniquevaluelist = []
    for i in list:
        if i.split()[index] in uniquevaluelist:
            nothing = 1
        else:
            uniquevaluelist.append(i.split()[index])
    return uniquevaluelist

def getdirs(path,tag,line):
    linelength = len(line)
    dirlist = os.listdir(path)

    changevar = []
    for i in dirlist:
        if tag in i:
            index = i.find(line + '_') + linelength + 1
            changevar.append(i[index:])
    return changevar

def mkposcars(file):
    input = readfile(file)
    

    if input[0].split()[0] == 'all':
        superstruct = getdirs('vaspruns','str','str')
    elif '-' in input[0].split():
        superstruct = []
        for i in range(int(input[0].split()[0]),int(input[0].split()[2])+1):
            superstruct.append(str(i))
    else:
        superstruct = input[0].split()


    parentatoms = int(input[6])
    element2 = 8 + parentatoms
    
    elementone = input[1].split()[0]
    elementtwo = input[element2].split()[0]

    totatoms = int(input[6])   #find total number of atoms in parent unit cell
    latparindex=9+totatoms     #find location of lattice parameter

    firstlatpar = float(input[2])
    secondlatpar = float(input[latparindex])

    if firstlatpar > secondlatpar:
        largelatpar=firstlatpar
        smalllatpar=secondlatpar
        n=1
        p=0
    else:
        largelatpar=secondlatpar
        smalllatpar=firstlatpar
        n=0
        p=1


    
    for i in superstruct:
        os.system('makestr.new struct_enum.out ' + i)

        if int(i) > 9999:
            poscarread = open('vasp.0' + i, 'r')    
        elif int(i) > 999:
            poscarread = open('vasp.00' + i, 'r')    
        elif int(i) > 99:
            poscarread = open('vasp.000' + i, 'r')
        elif int(i) > 9: 
            poscarread = open('vasp.0000' + i, 'r')
        else:
            poscarread = open('vasp.00000' + i, 'r')

        lines = poscarread.readlines()
        poscarread.close()
        atoms = lines[5].split()
        lines[1] = str(vegardslaw(atoms[n],atoms[p],largelatpar,smalllatpar)) + '\n'

        totalatoms = int(atoms[0]) + int(atoms[1])

        if int(i) > 9999:
            poscarwrite = open('vasp.0' + i, 'w')
        elif int(i) > 999:
            poscarwrite = open('vasp.00' + i, 'w')
        elif int(i) > 99:
            poscarwrite = open('vasp.000' + i, 'w')
        elif int(i) > 9:
            poscarwrite = open('vasp.0000' + i, 'w')
        else:
            poscarwrite = open('vasp.00000' + i, 'w')
    
        poscarwrite.writelines(lines)
        poscarwrite.close()


    indexone = 8+totatoms
    indextwo = input.index('kpoints XxXxX\n')

    writefile('vasp.' + elementone,input[1:indexone])
    
    writefile('vasp.' + elementtwo,input[indexone:indextwo])


def strtocartesian(basisvecs,latticevecs):
    cartesianbasisvecs = []

    for k in basisvecs:
        vecone = round(float(k.split()[0]) * float(latticevecs[0].split()[0]) + float(k.split()[1]) * float(latticevecs[1].split()[0]) + float(k.split()[2]) * float(latticevecs[2].split()[0]),8)
        vectwo = round(float(k.split()[0]) * float(latticevecs[0].split()[1]) + float(k.split()[1]) * float(latticevecs[1].split()[1]) + float(k.split()[2]) * float(latticevecs[2].split()[1]),8)
        vecthree = round(float(k.split()[0]) * float(latticevecs[0].split()[2]) + float(k.split()[1]) * float(latticevecs[1].split()[2]) + float(k.split()[2]) * float(latticevecs[2].split()[2]),8)
        cartesianbasisvecs.append([vecone,vectwo,vecthree])

    return cartesianbasisvecs

def strtocartesian2d(basisvecs,latticevecs):
    cartesianbasisvecs = []

    for k in basisvecs:
        vecone = round(float(k.split()[0]) * float(latticevecs[0].split()[0]) + float(k.split()[1]) * float(latticevecs[1].split()[0]),2)
        vectwo = round(float(k.split()[0]) * float(latticevecs[0].split()[1]) + float(k.split()[1]) * float(latticevecs[1].split()[1]),2)
        cartesianbasisvecs.append([vecone,vectwo])

    return cartesianbasisvecs



def xcombinations(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xcombinations(items[:i]+items[i+1:],n-1):
                yield [items[i]]+cc

def xuniqueCombinations(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xuniqueCombinations(items[i+1:],n-1):
                yield [items[i]]+cc
            
def xselections(items, n):
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for ss in xselections(items, n-1):
                yield [items[i]]+ss

def xpermutations(items):
    return xcombinations(items, len(items))


def cubicsymops():
    temp = []
    latveccombs = []
    tuples = []


    for z in xselections([1,-1],3):
        tuples.append(z)


    for x in xpermutations([1,0,0]):
        if x not in temp:
            temp.append(x)

    for y in xpermutations(temp):
        if y not in latveccombs:
            latveccombs.append(y)


    list = []

    for i in tuples:
        for j in latveccombs:
            vec = [[j[x][0]*i[x],j[x][1]*i[x],j[x][2]*i[x]] for x in range(0,3)]
            if vec not in list:
                list.append(vec)


    return list

def squaresymops():
    temp = []
    latveccombs = []
    tuples = []


    for z in xselections([1,-1],2):
        tuples.append(z)


    for x in xpermutations([1,0]):
        if x not in temp:
            temp.append(x)

    for y in xpermutations(temp):
        if y not in latveccombs:
            latveccombs.append(y)


    list = []

    for i in tuples:
        for j in latveccombs:
            vec = [[j[x][0]*i[x],j[x][1]*i[x]] for x in range(0,2)]
            if vec not in list:
                list.append(vec)


    return list



def matrixmultiply(vec,matrix):
    vecone = matrix[0][0]*vec[0] + matrix[0][1]*vec[1] + matrix[0][2]*vec[2]
    vectwo = matrix[1][0]*vec[0] + matrix[1][1]*vec[1] + matrix[1][2]*vec[2]
    vecthree = matrix[2][0]*vec[0] + matrix[2][1]*vec[1] + matrix[2][2]*vec[2]
    
    return [vecone,vectwo,vecthree]

def matrixmultiply2d(vec,matrix):
    vecone = round(matrix[0][0]*vec[0] + matrix[0][1]*vec[1])
    vectwo = round(matrix[1][0]*vec[0] + matrix[1][1]*vec[1])
    
    
    return [vecone,vectwo]



def mapback(latvecs,vec,basisvecs):
    vecs = []
    for i in range(-3,3):
        for j in range(-3,3):
            for k in range(-3,3):
                vecone = [i*float(latvecs[0].split()[0]),i * float(latvecs[0].split()[1]),i * float(latvecs[0].split()[2])]
                vectwo = [j*float(latvecs[1].split()[0]),j * float(latvecs[1].split()[1]),j * float(latvecs[1].split()[2])]
                vecthree = [k*float(latvecs[2].split()[0]),k * float(latvecs[2].split()[1]),k * float(latvecs[2].split()[2])]
                newvector = [vecone[0] + vectwo[0] + vecthree[0] + vec[0],vecone[1] + vectwo[1] + vecthree[1] + vec[1],vecone[2] + vectwo[2] + vecthree[2] + vec[2]]
                vecs.append(newvector)

                
    index = 0
    returnvar = '????????'
    for l in basisvecs:
        for m in vecs:
            if l == m:
                print basisvecs
                print 'basisvecs'
                print index + 1
                print 'INDEX'
                returnvar = basisvecs[index + 1]
                print returnvar
        index = index + 1

    return returnvar

def mapback2d(latvecs,vec,basisvecs,upperlimit):
    vecs = []
    for j in range(-3,2):
        for k in range(-3 ,2):
            vecone = [j * float(latvecs[0].split()[0]),j * float(latvecs[0].split()[1])]
            vectwo = [k * float(latvecs[1].split()[0]),k * float(latvecs[1].split()[1])]
             
            newvector = [vecone[0] + vectwo[0] + vec[0],vecone[1] + vectwo[1] + vec[1]]
            vecs.append(newvector)

                
    index = 0
    returnvar = '????????'
    for l in basisvecs:
        for m in vecs:
            if l == m:
                returnvar = m
        index = index + 1

    return returnvar




def twovectoradd(vectorone,vectortwo):
    return [vectorone[0] + vectortwo[0],vectorone[1] + vectortwo[1],vectorone[2] + vectortwo[2]]

def twovectoradd2d(vectorone,vectortwo):
    return [vectorone[0] + vectortwo[0],vectorone[1] + vectortwo[1]]


def threevectoradd(vectorone,vectortwo,vectorthree):
    return [vectorone[0] + vectortwo[0] + vectorthree[0],vectorone[1] + vectortwo[1] + vectorthree[1],vectorone[2] + vectortwo[2] + vectorthree[2]]

def scalartimesvector(scalar,vector):
    return [scalar * float(vector.split()[0]),scalar * float(vector.split()[1]),scalar * float(vector.split()[2])]
                    
def scalartimesvector2d(scalar,vector):
    return [scalar * float(vector.split()[0]),scalar * float(vector.split()[1])]


def eliminateduplicates2d(figurelist,uniquefigs,parent,symops,figtype):
    for j in figurelist:
            test = 'n'
            if j not in uniquefigs:
                for k in parent:
                    for l in symops:
                        trans = [twovectoradd2d(j[x],k) for x in range(0,figtype)]
                        rot = [matrixmultiply2d(trans[x],l) for x in range(0,figtype)]
                        for f in uniquefigs:
                            if figtype == 2 and rot[0] in f and rot[1] in f:
                                test = 'y'
                            if figtype == 3 and rot[0] in f and rot[1] in f and rot[2] in f:
                                test = 'y'
                            if figtype == 4 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f:
                                test = 'y'
                            if figtype == 5 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f:
                                test = 'y'
                            if figtype == 6 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and rot[5] in f:
                                test = 'y'
                            if figtype == 7 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and rot[5] in f and rot[6] in f:
                                test = 'y'
                            if figtype == 8 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and rot[5] in f and rot[6] in f and rot[7] in f:
                                test = 'y'
                                
                if test == 'n':
                    uniquefigs.append(j)

    return uniquefigs

def eliminateduplicates(figurelist,uniquefigs,parent,symops,figtype):
    for j in figurelist:
            test = 'n'
            if j not in uniquefigs:
                for k in parent:
                    for l in symops:
                        trans = [twovectoradd(j[x],k) for x in range(0,figtype)]
                        rot = [matrixmultiply(trans[x],l) for x in range(0,figtype)]
                        for f in uniquefigs:
                            if figtype == 2 and rot[0] in f and rot[1] in f:
                                test = 'y'
                            if figtype == 3 and rot[0] in f and rot[1] in f and rot[2] in f:
                                test = 'y'
                            if figtype == 4 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f:
                                test = 'y'
                            if figtype == 5 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f:
                                test = 'y'
                            if figtype == 6 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and rot[5] in f:
                                test = 'y'
                            if figtype == 7 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and rot[5] in f and rot[6] in f:
                                test = 'y'
                            if figtype == 8 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and rot[5] in f and rot[6] in f and rot[7] in f:
                                test = 'y'

                                
                if test == 'n':
                    uniquefigs.append(j)

    return uniquefigs


def findfigures2d(figurelist,structurefigs,parent,symops,figtype):
    templist = []
    for j in structurefigs:
        for k in parent:
            for l in symops:
                trans = [twovectoradd2d(j[x],k) for x in range(0,figtype)]
                rot = [matrixmultiply2d(trans[x],l) for x in range(0,figtype)]
                for f in figurelist:
                    if figtype == 2 and rot[0] in f and rot[1] in f and f not in templist:
                        templist.append(f)
                    if figtype == 3 and rot[0] in f and rot[1] in f and rot[2] in f and f not in templist:
                        templist.append(f)
                    if figtype == 4 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and f not in templist:
                        templist.append(f)
                    if figtype == 5 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and f not in templist:
                        templist.append(f)                        
                    if figtype == 6 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and rot[5] in f and f not in templist:
                        templist.append(f)                        
                    if figtype == 7 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and rot[5] in f and rot[6] in f and f not in templist:
                        templist.append(f)                        
                    if figtype == 8 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and rot[5] in f and rot[6] in f and rot[7] in f and f not in templist:
                        templist.append(f)                        
                                
                    

    return templist

def findfigures(figurelist,structurefigs,parent,symops,figtype):
    templist = []
    for j in structurefigs:
        for k in parent:
            for l in symops:
                trans = [twovectoradd(j[x],k) for x in range(0,figtype)]
                rot = [matrixmultiply(trans[x],l) for x in range(0,figtype)]
                for f in figurelist:
                    if figtype == 2 and rot[0] in f and rot[1] in f and f not in templist:
                        templist.append(f)
                    if figtype == 3 and rot[0] in f and rot[1] in f and rot[2] in f and f not in templist:
                        templist.append(f)
                    if figtype == 4 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and f not in templist:
                        templist.append(f)
                    if figtype == 5 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and f not in templist:
                        templist.append(f)                        
                    if figtype == 6 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and rot[5] in f and f not in templist:
                        templist.append(f)                        
                    if figtype == 7 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and rot[5] in f and rot[6] in f and f not in templist:
                        templist.append(f)                        
                    if figtype == 8 and rot[0] in f and rot[1] in f and rot[2] in f and rot[3] in f and rot[4] in f and rot[5] in f and rot[6] in f and rot[7] in f and f not in templist:
                        templist.append(f)                        
                                
                    

    return templist


def surfacedistancecheck(l):
    if int(l) > 9999:
        poscar = readfile('vasp.0' + l)
    elif int(l) > 999:
        poscar = readfile('vasp.00' + l)
    elif int(l) > 99:
        poscar = readfile('vasp.000' + l)
    elif int(l) > 9:
        poscar = readfile('vasp.0000' + l)
    else:
        poscar = readfile('vasp.00000' + l)
    lvtwo = poscar[3].split()[1] + ' ' + poscar[3].split()[2]
    lvthree = poscar[4].split()[1] + ' ' + poscar[4].split()[2]
    lvecs = [lvtwo,lvthree]
    numofbasis = int(poscar[5].split()[0])
    basisvecs = []
    for k in range(7,7+numofbasis):
        basisvecs.append(poscar[k].split()[1] + ' ' + poscar[k].split()[2])
    cartesianvecs=strtocartesian2d(basisvecs,lvecs)
    lengths = []
    for m in cartesianvecs:
        for l in cartesianvecs:
            diff = [float(m[0])-float(l[0]),float(m[1])-float(l[1])]
            distance = magnitude2d(diff)
            if distance not in lengths and l != m:
                lengths.append(distance)
    for r in lengths:
        if float(r) <= 0.5:
            status = 'dont make'
            break
        else:
            status = 'make it so'
    if len(lengths) == 0:
        status = 'make it so'

    return status
