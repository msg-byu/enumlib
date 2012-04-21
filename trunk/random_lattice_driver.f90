PROGRAM Reducedriver
use num_types
use vector_matrix_utilities
use rational_mathematics
use cwrapper
implicit none

real(dp) rand(9), R(3,3), Out(3,3)
integer L(3,3), H(3,3),S(3,3), T1(2,2), T2(2,2), T3(2,2)
integer i,j,sz

sz=10
call random_seed()
!do i=1,10
   call random_number(rand)
   print *, "starting new matrix"
   R = reshape((rand*sz)-sz/2,(/3,3/))
   !call reduce_to_shortest_basis(R,out,1e-12_dp)
   call aflow_reduce_to_shortest_basis(R,out,1e-12_dp)
   !H = reshape((/(/1,0,0/),(/0,3,0/),(/0,-2,2/)/),(/3,3/))
   !call SmithNormalForm(H,L,S,R)
!enddo

do i = 1,3
   write(*,'(2(3(f7.3,1x)))') R(:,i), Out(:,i)
enddo


ENDPROGRAM Reducedriver
