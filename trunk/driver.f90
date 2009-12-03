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

character(80) title
logical fullLab
integer, pointer :: label(:,:)
integer, pointer :: digit(:)
integer, pointer :: equivalencies(:)

call read_input(title,LatDim,parLV,nD,d,k,equivalencies,nMin,nMax,eps,fullLab,label,digit) ! Read in parent lattice vectors, etc.
if (LatDim==3) then; latTyp='b';else;latTyp='s';endif
! With test case 006 there is a problem with the original code. label rotation fails...
! call generate_derivative_structures(title, parLV,nD,d,k,nMin,nMax,latTyp,eps,fullLab)

!!allocate(label(3,4),digit(4))
!!label(:,1) = (/3,1,0/)
!!label(:,2) = (/0,2,0/)
!!label(:,3) = (/1,4,0/)
!!label(:,4) = (/0,1,2/)
!!digit =(/2,2,2,3/)
!!
!!call mixed_radix_counter(label,digit)
!!
call gen_multilattice_derivatives(title, parLV,nD,d,k,nMin,nMax,latTyp,eps,fullLab,&
         label,digit,equivalencies,conc_check=.false.)

END PROGRAM driver
