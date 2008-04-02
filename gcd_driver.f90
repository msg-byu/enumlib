program gcd_driver
use vector_matrix_utilities
use rational_mathematics
implicit none

real x(3,3)
integer c
integer, dimension(3,3) ::  H, M, A, B
integer i,j,i1,i2,i3
!read(*,*) a,b
!print *, "here"
!write(*,*) "gcd of ",a,b
!write(*,*) gcd(a,b)
!print *,"there"
!
!read(*,*) x
!write(*,*) "gcd of",x
!write(*,*) gcd(x)

H = reshape((/2,0,3, 0,5,0, 4,6,2/),(/3,3/))

call system_clock(i1,i2,i3)
call random_seed(put=(/i1/))
do
   call random_number(x)
   H = nint(x*10)-5
   if (determinant(H)/=0) exit
enddo
!H = reshape((/-3,-5,-1, -4,-3,-3, 0,1,3/),(/3,3/))
!H = reshape((/-2,-2,-2, 0,1,2, 3,2,3/),(/3,3/))
!H = reshape((/-4,2,4, 0,-4,-3, 4,-4,3/),(/3,3/))
!H = reshape((/-2,-2,-2, -4,-4,4, -2,-1,-4/),(/3,3/))
!H = reshape((/5,3,-3, -4,4,4, 1,-4,-1/),(/3,3/))

print *,"Determinant is",determinant(H)
call SmithNormalForm(H,A,M,B)
print *,"Determinant is",determinant(H)

M = matmul(matmul(A,H),B)
do i = 1,3
   write(*,'(3(3x,3i3))') M(i,:)
enddo

endprogram gcd_driver
