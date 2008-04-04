PROGRAM driver
use num_types 
use derivative_structure_generator
use io_utilities
implicit none

integer nMin, nMax ! Numbers of various things
integer k ! Number of colors/label types (i.e., binary, ternary, etc.)
integer LatDim ! 2D or 3D parent lattice?
real(dp)  eps
real(dp) :: parLV(3,3)

character(80) title

call read_input(title,LatDim,parLV,k,nMin,nMax,eps) ! Read in parent lattice vectors, etc.
 
call generate_derivative_structures(title, parLV,k,nMin,nMax,LatDim,eps)

END PROGRAM driver
