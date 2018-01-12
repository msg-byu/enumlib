LOG		= compile.log
#FFLAGS		= -c
.SILENT:

# sources
WSRCF90		=  \
		f90wrap_HNF_profiler.f90 \

SRCF90		=  \
		../aux_src/HNF_profiler.f90 \

OBJSF90		= $(SRCF90:.f90=.o)

LBD1 = ../symlib/src
LBD2 = ../src
LIBS =  $(LBD2)/libenum.a $(LBD1)/libcomparestructs.a $(LBD1)/libutils.a \
         $(LBD1)/librational.a $(LBD1)/libcombinatorics.a $(LBD1)/libsym.a 

#IMPORTANT: remember to include the *.o and libraries in the call to create
#the shared library!
all: pre_comp $(OBJSF90)
	-rm _supercell.so
	echo "Creating shared library _supercell.so"
	-f2py-f90wrap -c -m _supercell $(WSRCF90) *.o ${LIBS} -I${LDB1} -I${LDB2} -I.>> /dev/null
	echo "done."

pre_comp:
	+$(MAKE) -C ../symlib/src
	+$(MAKE) -C ../src/

$(OBJSF90): %.o: %.f90
	echo "Compiling: $^"
	-$(F90) -fPIC -I${LBD1} -I${LBD2} $^ -o $@ ${LIBS}>> $(LOG) 2>> $(LOG)
	-$(F90) -fPIC -c -I${LBD1} -I${LBD2} $^
	echo "done."

clean:
	+$(MAKE) clean -C ../symlib/src
	+$(MAKE) clean -C ../src
	-rm *.o *.mod *.so compile.log
