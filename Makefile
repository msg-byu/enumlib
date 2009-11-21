# Set up conditional compilation so that one Makefile can handle
# multiple platforms/compilers. Branch according to an environmental
# variable F90. I wish someone would show me a better way of doing this.
#
LBDR = ../../celib/trunk
FOUND = false
ifeq (${F90},gfortran)  # gfortran compiler
  FFLAGS =  -g -fbounds-check -warn  -I${LBDR}
  FOUND = true
endif
ifeq (${F90},ifc)  # Intel compiler
  FFLAGS =  -g -error-limit 7 -traceback -check bounds -warn  -I${LBDR}
  FOUND = true
endif

ifeq (${F90},ifort)  # Intel compiler
  ifeq (${DEBUG},false)
     FFLAGS =  -O3 -xW -I${LBDR} 
     FOUND = true
  endif
  ifeq (${DEBUG},true)
#  F90 =  /opt/intel/fc/10.0.016/bin/ifort
     FFLAGS =  -g -debug -error-limit 7 -heap-arrays -traceback -check bounds -warn -e95 -I${LBDR} 
#-prof-use -prof-dir .
     FOUND = true
  endif
  ifeq (${DEBUG},)
     FFLAGS =  -g -debug -error-limit 7 -heap-arrays -traceback -check bounds -warn -I${LBDR}  
     FOUND = true
  endif

endif


ifeq (${F90},xlf90) # IBM compiler
  FFLAGS = -g -C -qsuffix=f=f90  -I${LBDR}
  FOUND = true
endif
ifeq (${F90},f95) # Absoft PPC compiler
#  FFLAGS =  -profile -p ${LBDR} 
  FFLAGS = -g -Rb -Rc  -nodefaultmod -p ${LBDR} #-ea
# B80  show entry in subprograms ; Rb bounds; Rc array conformance;
# z2 warning level
# -O3 optimization
# -ea stop after one error 
  FOUND = true
endif
ifeq (${F90},)  # If the variable isn't set, make sure compilation bombs
error:
	echo Error environmental variable "F90" not set!
endif
ifneq (${FOUND},true) # If it's set but no matching flags then bomb
error:	
	echo Error: makefile doesn\'t have flags for this compiler
endif

SRC = sorting.f90 enumeration_types.f90 io_utils.f90 labeling_related.f90 \
      enumeration_utilities.f90 derivative_structure_generator.f90 

OBJS = ${SRC:.f90=.o}
LIBS =  ${LBDR}/libcomparestructs.a ${LBDR}/libutils.a ${LBDR}/libsym.a \
         ${LBDR}/librational.a ${LBDR}/libcombinatorics.a 

.SUFFIXES :  
.SUFFIXES :  .f .f90 .f95 .o



libenum.a: ${OBJS}
	ar ru $@ $?
	ranlib  $@

all: libenum.a enum.x compare.x

multienum.x: ${OBJS} driver.o
	${F90} ${LDFLAGS} -o $@ ${OBJS} driver.o ${LIBS}

compare.x: compare.o
	${F90} ${LDFLAGS} -o $@ compare.o ${LIBS}

2Dplot.x: make2Dplot.o splot.o
	${F90} ${LDFLAGS} -o $@ splot.o make2Dplot.o ${LIBS}

makestr.x: makeStr.o
	${F90} ${LDFLAGS} -o $@ makeStr.o ${LIBS}

makeperovstr.x: makePerovStr.o
	${F90} ${LDFLAGS} -o $@ makePerovStr.o ${LIBS}

makestructin.x: makeStrIn.o
	${F90} ${LDFLAGS} -o $@ makeStrIn.o ${LIBS}

.f95.o : 
	${F90} ${FFLAGS} -c $<
.f90.o : 
	${F90} ${FFLAGS} -c $<
.f.o : 
	${F90} -c $<



CLEAN  = *.o *.mod *.a
clean : 
	rm -f ${CLEAN}
clobber : 
	rm -f  ${CLEAN}  *~ \#*
	make
