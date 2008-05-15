# Set up conditional compilation so that one Makefile can handle
# multiple platforms/compilers. Branch according to an environmental
# variable F90. I wish someone would show me a better way of doing this.
#
ifeq (${F90},ifc)  # Intel compiler
  LBDR = ../celib
  FFLAGS =  -g -error-limit 7 -traceback -check bounds -warn  -I${LBDR}
  FOUND = true
endif
ifeq (${F90},ifort)  # Intel compiler
  F90 = /opt/intel/fc/10.0.016/bin/ifort
  LBDR = ../celib
  FFLAGS =  -g -debug -error-limit 7 -traceback -check bounds -warn -e95  -I${LBDR}
  FOUND = true
endif
ifeq (${F90},xlf90) # IBM compiler
  LBDR = ../celib
  FFLAGS = -g -C -qsuffix=f=f90  -I${LBDR}
  FOUND = true
endif
ifeq (${F90},f95) # Absoft PPC compiler
  LBDR = ../celib
  FFLAGS = -g -Rb -Rc -z2 -et -nodefaultmod -p ${LBDR} 
#  FFLAGS = -O3 -nodefaultmod -p ../guslib/ #-ea
# B80  show entry in subprograms ; Rb bounds; Rc array conformance;
# z2 warning level
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

SRC = enumeration_types.f90 labeling_related.f90 derivative_structure_generator.f90 \
	io_utils.f90 
OBJS = ${SRC:.f90=.o}
LIBS =  ${LBDR}/libcomparestructs.a ${LBDR}/libutils.a ${LBDR}/libsym.a \
         ${LBDR}/librational.a ${LBDR}/libcombinatorics.a 

.SUFFIXES :  
.SUFFIXES :  .f .f90 .f95 .o



libenum.a: ${OBJS}
	ar ru $@ $?
	ranlib  $@

all: libenum.a enum.x compare.x

enum.x: ${OBJS} driver.o
	${F90} ${LDFLAGS} -o $@ ${OBJS} driver.o ${LIBS}

compare.x: compare.o
	${F90} ${LDFLAGS} -o $@ compare.o ${LIBS}


.f95.o : 
	${F90} ${FFLAGS} -c $<

.f90.o : 
	${F90} ${FFLAGS} -c $<

CLEAN  = *.o *.mod *.a
clean : 
	rm -f ${CLEAN}
	make
clobber : 
	rm -f  ${CLEAN}  *~ \#*
	make