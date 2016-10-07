PROGRAM driver
use num_types 
use derivative_structure_generator
use io_utils
implicit none

integer nMin, nMax ! Numbers of various things
integer k ! Number of colors/label types (i.e., binary, ternary, etc.)
integer LatDim ! 2D or 3D parent lattice?
integer nD ! Number of sites in the basis (i.e., number of points in the multilattice)
real(dp), pointer :: d(:,:) => null()
character(1) :: latTyp
real(dp)  eps
real(dp) :: parLV(3,3)

character(80) title, fname
logical fullLab,concCheck
integer, pointer :: cRange(:,:)
integer, pointer :: label(:,:)
integer, pointer :: digit(:)
integer, pointer :: equivalencies(:)

if (iargc()>=1) then
   call getarg(1,fname)
else
   fname = "struct_enum.in"
endif
call read_input(title,LatDim,parLV,nD,d,k,equivalencies,nMin,nMax,eps&
     &,fullLab,label,digit,fname,cRange,concCheck) ! Read in parent lattice vectors, etc.
if (LatDim==3) then; latTyp='b';else;latTyp='s';endif
call gen_multilattice_derivatives(title, parLV,nD,d,k,nMin,nMax,latTyp,eps,fullLab,&
         label,digit,equivalencies,concCheck,cRange, polya=.true.)

END PROGRAM driver
