LCLB = ./
LBDR = ../
FOUND = false
F90=gfortran

ifeq (${F90},gfortran)  # gfortran compiler
	FFLAGS =  -g -fbounds-check -Wall -ffree-line-length-none -I${LCLB} -I${LBDR}
	FOUND = true
endif

ifeq (${F90},ifc)  # Intel compiler
	FFLAGS =  -g -error-limit 7 -traceback -check bounds -warn  -I${LCLB} -I${LBDR}
	FOUND = true
endif

ifeq (${F90},ifort)  # Intel compiler
	ifeq (${DEBUG},false)
		FFLAGS =  -O3 -I${LBDR}
		FOUND = true
	endif
 	ifeq (${DEBUG},true)
 		FFLAGS =  -g -debug -error-limit 7 -heap-arrays -traceback -check bounds -warn -I${LBDR}
		FOUND = true
	endif
endif

SRC = derivative_structure_gen_test.f90 
OBJ = ${SRC:.f90=.o}\

all: bin

bin : $(OBJ) 
	${F90} test_driver.f90 -o test.x $(OBJ)
	chmod ugo+x test.x

derivative_structure_gen_test.o : derivative_structure_gen_test.f90
	${F90} -c derivative_structure_gen_test.f90 -I ../
	
	
	
clean :
	- echo "Removing generated files"
	- rm -rf *.o
	- rm -rf *.mod
	- rm -rf *.x