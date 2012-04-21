!***************************************************************************************************
! This wrapper provides an easy interface for a C++ program to call something in the 
MODULE cwrapper
use num_types
use numerical_utilities
use vector_matrix_utilities
use enumeration_utilities

implicit none
private
public aflow_reduce_to_shortest_basis

CONTAINS

!***************************************************************************************************
! Make a call to the generator for derivative structures
SUBROUTINE aflow_reduce_to_shortest_basis(A,B,eps)
real(dp) :: A(3,3), B(3,3), eps
print*,"eps",eps
if (equal(determinant(A),0._dp,1.e-12_dp)) stop "Matrix is singular in aflow_reduce_to_shortest_basis"
call reduce_to_shortest_basis(A,B,eps)
END SUBROUTINE aflow_reduce_to_shortest_basis

SUBROUTINE aflow_reduce_to_shortest_basis2(rA,cA,A,rB,cB,B,eps)
real(dp) :: A(rA,cA), eps
real(dp), pointer :: B(:,:)
!real(dp), intent(out) ::  B(:,:)
integer  :: rA, cA, rB, cB, i
print*,"eps",eps
allocate(B(3,3))
write(*,'("A",/,3(f7.4,1x))') (A(:,i),i=1,3)
B = reshape((/1,2,3,4,5,6,7,8,9/),(/3,3/))
write(*,'("B",/,3(f7.4,1x))') (B(:,i),i=1,3)
!write(*,'("B",/,3(f7.4,1x))') (B(:,i),i=1,3)
if (equal(determinant(A),0._dp,1.e-12_dp)) stop "Matrix is singular in aflow_reduce_to_shortest_basi&
     &s2"
!call reduce_to_shortest_basis(A,B,eps)
rB=size(B,1); cB=size(B,2)
END SUBROUTINE aflow_reduce_to_shortest_basis2

END MODULE cwrapper

