program snf_driver
use vector_matrix_utilities
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
!print *,i1,i2,i3
call random_seed(put=(/i1/))
do
   call random_number(x)
   H = nint(x*10)-5
   if (determinant(H)/=0) exit
enddo
print *,"Determinant is",determinant(H)
call smith_normal_form(H,3,3,3,A,3,B,3)
do i = 1,3
   write(*,'(3(3x,3i3))') A(i,:),M(i,:),B(i,:)
enddo
print *,"Determinant is",determinant(H)


endprogram snf_driver
