# Set up conditional compilation so that one Makefile can handle
# multiple platforms/compilers. Branch according to an environmental
# variable F90. I wish someone would show me a better way of doing this.
ifeq (${F90},ifc)  # Intel compiler
  LBDR = ../celib
  FFLAGS =  -g -error-limit 7 -traceback -check bounds -warn  -I${LBDR}
  FOUND = true
endif
ifeq (${F90},ifort)  # Intel compiler
  LBDR = ../celib
  FFLAGS =  -g -error-limit 7 -traceback -check bounds -warn -e95  -I${LBDR}
  FOUND = true
endif
ifeq (${F90},xlf90) # IBM compiler
  LBDR = ../celib
  FFLAGS = -g -C -qsuffix=f=f90  -I${LBDR}
  FOUND = true
endif
ifeq (${F90},f95) # Absoft PPC compiler
  LBDR = ../celib
  FFLAGS = -g -Rb -Rc  -z2 -nodefaultmod -p ${LBDR} 
#  FFLAGS = -O3 -nodefaultmod -p ../guslib/ #-ea
# B80  show entry in subprograms ; Rb bounds; Rc array conformance; 
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

# Main parts of the makefile are here
EXENAME = enum.x
SRC = io_utilities.f90 crystal_types.f90 labeling_related.f90 \
      derivative_structure_generator.f90 driver.f90 
OBJS = ${SRC:.f90=.o}
LIBS =  ${LBDR}/libcomparestructs.a ${LBDR}/libutils.a ${LBDR}/libsym.a \
         ${LBDR}/librational.a ${LBDR}/libcombinatorics.a 

# Delete all default suffix rules
.SUFFIXES :  
.SUFFIXES :  .f90 .f .o .mod

${EXENAME}: ${OBJS}    
	${F90} ${LDFLAGS} -o ${EXENAME} ${OBJS} ${LIBS}

.f90.o : 
	${F90} ${FFLAGS} -c $<
.f.o : 
	${F90} -c $<

poly_test: ${OBJS} poly_driver.o
	${F90} ${LDFLAGS} -o ${EXENAME} ${OBJS} poly_driver.o ${LIBS}

makestr: makeStr.o
	${F90} ${LDFLAGS} -o makestr.x makeStr.o ${LIBS}

2Dplot: make2Dplot.o splot.o
	${F90} ${LDFLAGS} -o 2D.x splot.o make2Dplot.o ${LIBS}
clean: 
	rm -f *.o ${EXENAME} *.mod
	make
clobber: 
	rm -f ${OBJS} ${EXENAME} *.mod *~ \#*
