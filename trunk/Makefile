# Set up conditional compilation so that one Makefile can handle
# multiple platforms/compilers. Branch according to an environmental
# variable F90. I wish someone would show me a better way of doing this.
#
LBDR = ../../celib/trunk
FOUND = false
ifeq (${F90},gfortran)  # gfortran compiler
  FFLAGS = -fPIC -g -fbounds-check -Wall -ffree-line-length-none -I${LBDR} 
  FOUND = true
endif

ifeq (${F90},ifc)  # Intel compiler
  FFLAGS = -fPIC -g -error-limit 7 -traceback -check bounds -warn  -I${LBDR}
  FOUND = true
endif

ifeq (${F90},ifort)  # Intel compiler
  ifeq (${DEBUG},false)
     FFLAGS =  -fPIC -O3 -I${LBDR} 
     FOUND = true
  endif
  ifeq (${DEBUG},true)
#  F90 =  /opt/intel/fc/10.0.016/bin/ifort
     FFLAGS =  -fPIC -g -debug -error-limit 7 -heap-arrays -traceback -check bounds -warn -I${LBDR} 
#-prof-use -prof-dir .
     FOUND = true
  endif
  ifeq (${DEBUG},)
     FFLAGS =  -g -fPIC -debug -error-limit 7 -heap-arrays -traceback -check bounds -warn -I${LBDR}  
     FOUND = true
  endif

endif


ifeq (${F90},xlf90) # IBM compiler
  FFLAGS = -g -C -fPIC -qsuffix=f=f90  -I${LBDR}
  FOUND = true
endif
ifeq (${F90},f95) # Absoft PPC compiler
#  FFLAGS =  -profile -p ${LBDR} 
  FFLAGS = -g -Rb -Rc  -fPIC -nodefaultmod -p ${LBDR} #-ea
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
      derivative_structure_generator.f90 enumeration_utilities.f90

OBJS = ${SRC:.f90=.o}
LIBS =  ${LBDR}/libcomparestructs.a ${LBDR}/libutils.a ${LBDR}/libsym.a \
         ${LBDR}/librational.a ${LBDR}/libcombinatorics.a 

.SUFFIXES :  
.SUFFIXES :  .f .f90 .f95 .o



libenum.a: ${OBJS}
	ar ru $@ $?
	ranlib  $@

all: libenum.a multienum.x find_structure_in_list.x 2Dplot.x \
     compare_two_enum_files.x  
multienum.x: ${OBJS} driver.o
	${F90} ${LDFLAGS} -o $@ ${OBJS} driver.o ${LIBS}

2Dplot.x: make2Dplot.o splot.o
	${F90} ${LDFLAGS} -o $@ splot.o make2Dplot.o ${LIBS}

makestr.x: makeStr.o
	${F90} ${LDFLAGS} -o $@ makeStr.o libenum.a ${LIBS} 

find_structure_in_list.x: ${OBJS} find_structure_in_list.o
	${F90} ${LDFLAGS} -o $@ enumeration_utilities.o find_structure_in_list.o libenum.a ${LIBS}
compare_two_enum_files.x: ${OBJS} compare_two_enum_files.o
	${F90} ${LDFLAGS} -o $@ ${OBJS} compare_two_enum_files.o libenum.a ${LIBS}
convert_structures_to_enumformat.x: ${OBJS} convert_structures_to_enumformat.o
	${F90} ${LDFLAGS} -o $@  enumeration_utilities.o convert_structures_to_enumformat.o libenum.a ${LIBS}

makestr.2d: makeStr2d.o
	${F90} ${LDFLAGS} -o $@ makeStr2d.o libenum.a ${LIBS} 

makeperovstr.x: makePerovStr.o
	${F90} ${LDFLAGS} -o $@ makePerovStr.o ${LIBS}

makestructin.x: makeStrIn.o
	${F90} ${LDFLAGS} -o $@ makeStrIn.o libenum.a ${LIBS}

randReduceTest.x: random_lattice_driver.o
	${F90} ${LDFLAGS} -o $@ random_lattice_driver.o libenum.a ${LIBS}


.f95.o : 
	${F90} ${FFLAGS} -c $<
.f90.o : 
	${F90} ${FFLAGS} -c $<
.f.o : 
	${F90} -c $<



CLEAN  = *.o *.mod *.a *.x
clean : 
	rm -f ${CLEAN}
clobber : 
	rm -f  ${CLEAN}  *~ \#*
	make
