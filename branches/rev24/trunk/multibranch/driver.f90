PROGRAM driver
use num_types 
use derivative_structure_generator
use io_utils
implicit none

integer nMin, nMax ! Numbers of various things
integer k ! Number of colors/label types (i.e., binary, ternary, etc.)
integer LatDim ! 2D or 3D parent lattice?
character(1) :: latTyp
real(dp)  eps
real(dp) :: parLV(3,3)

character(80) title
logical fullLab

call read_input(title,LatDim,parLV,k,nMin,nMax,eps,fullLab) ! Read in parent lattice vectors, etc.
if (LatDim==3) then; latTyp='b';else;latTyp='s';endif
call generate_derivative_structures(title, parLV,k,nMin,nMax,latTyp,eps,fullLab)

END PROGRAM driver
