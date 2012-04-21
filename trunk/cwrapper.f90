!***************************************************************************************************
! This wrapper provides an easy interface for a C++ program to call something in the 
MODULE cwrapper
use num_types
use numerical_utilities
use vector_matrix_utilities
use enumeration_utilities

implicit none
private
public Ccall2_reduce

CONTAINS

!***************************************************************************************************
! Make a call to the generator for derivative structures
SUBROUTINE aflow_reduce_to_shortest_basis(A,B,eps)
real(dp) :: A(3,3), B(3,3), eps
if (equal(determinant(A),0._dp,1.e-12_dp)) stop "Matrix is singular in Ccall2_reduce"
call reduce_to_shortest_basis(A,B,eps)
END SUBROUTINE aflow_reduce_to_shortest_basis

END MODULE cwrapper

