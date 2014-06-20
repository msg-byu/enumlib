PROGRAM deriv_driver
use num_types
use derivative_structures
use io_utils

type(enum_list) :: list

integer nMin, nMax ! Numbers of various things
integer k ! Number of colors/label types (i.e., binary, ternary, etc.)
integer LatDim ! 2D or 3D parent lattice?
integer nD ! Number of sites in the basis (i.e., number of points in the multilattice)
real(dp), pointer :: d(:,:) => null()
character(1) :: latTyp
real(dp)  eps
real(dp) :: parLV(3,3)

character(len=:), allocatable :: title, fname
character(5000) :: buffer
logical fullLab,concCheck
integer, pointer :: cRange(:,:)
integer, pointer :: label(:,:)
integer, pointer :: digit(:)
integer, pointer :: equivalencies(:)


if (iargc()>=1) then
   call getarg(1,buffer)
   fname = trim(buffer)
else
   fname = "struct_enum.out"
endif

call list%read_in_list(fname)


END PROGRAM deriv_driver
