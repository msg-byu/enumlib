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

call read_input(title,LatDim,parLV,nD,d,k,nMin,nMax,eps,fullLab) ! Read in parent lattice vectors, etc.
if (LatDim==3) then; latTyp='b';else;latTyp='s';endif
! With test case 006 there is a problem with the original code. label rotation fails...
! call generate_derivative_structures(title, parLV,nD,d,k,nMin,nMax,latTyp,eps,fullLab)

call gen_multilattice_derivatives(title, parLV,nD,d,k,nMin,nMax,latTyp,eps,fullLab,conc_check=.false.)

END PROGRAM driver
