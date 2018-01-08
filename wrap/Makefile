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
LIBS =  $(LBD1)/libcomparestructs.a $(LBD1)/libutils.a $(LBD1)/libsym.a \
         $(LBD1)/librational.a $(LBD1)/libcombinatorics.a $(LBD2)/libenum.a

#IMPORTANT: remember to include the *.o and libraries in the call to create
#the shared library!
all: $(OBJSF90)
	-rm _supercell.so
	echo "Creating shared library _supercell.so"
	-f2py-f90wrap -c -m _supercell $(WSRCF90) *.o ${LIBS} -I. >> /dev/null
	echo "done."

$(OBJSF90): %.o: %.f90
	echo "Compiling: $^"
	-$(F90) -fPIC -I${LBD1} -I${LBD2} $^ -o $@ ${LIBS}>> $(LOG) 2>> $(LOG)
	-$(F90) -fPIC -c -I${LBD1} -I${LBD2} $^
	echo "done."

clean:
	-rm *.o *.mod *.so
