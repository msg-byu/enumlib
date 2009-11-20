! This program reads in a file in struct_enum.out format and writes out
! a VASP-like output for each structure in the list, in a single file, generating a structures.in
! file that uncle can read. This is a good way to capture the correlations for a list of structures.
PROGRAM makeStrIn
use num_types
use vector_matrix_utilities
implicit none
character(80) fname, title, labeling
integer ioerr,  z1, z2, z3, ic, i, ilab, pgOps, nD, status
integer k, strN, sizeN, nAt, diag(3), a,b,c,d,e,f, HNF(3,3), L(3,3), hnfN
real(dp) :: p(3,3), sLV(3,3), Sinv(3,3), Binv(3,3), sLVorig(3,3), v(3), sLVinv(3,3)
real(dp) :: eps
real(dp), allocatable :: aBas(:,:), dvec(:,:)
character(1) bulksurf

character(40) dummy 
read(*,*) fname
open(11,file=fname,status='old',iostat=ioerr)
if(ioerr/=0)then; write(*,'("Input file doesn''t exist:",a80)') trim(fname);endif
! Read in the title from the struct_enum.out file
read(11,'(a80)') title; title = trim(title)

! Read in surf/bulk mode marker
read(11,'(a1)') bulksurf

! Read in the parent lattice vectors
do i = 1,3; read(11,*) p(:,i); enddo
call matrix_inverse(p,Binv)
!call matrix_inverse(p,pinv,err); if (err) stop "Parent lattice vectors are coplanar (makeStrIn)"

! Read in the number of d-vectors, then read each one in
read(11,*) nD
allocate(dvec(3,nD))
do i = 1,nD; read(11,*) dvec(:,i); enddo

! Read in the number of labels, i.e., binary, ternary, etc.
read(11,'(i2)') k
read(11,*) ! Skip one line
read(11,*) eps
i = 0
do 
   i = i + 1
   read(11,*) dummy
   if(dummy=="start") exit
   if (i>20) stop "Didn't find the ""start"" tag in the input file"
enddo 

open(12,file="structures.in") ! The output goes in this file (despite the "in" name)
write(12,'("peratom")')
write(12,'("noweights")')
do ! loop until end of file
   ! Read in the info for the given structure
   read(11,*,iostat=status) strN, hnfN,sizeN, nAt, pgOps, diag, a,b,c,d,e,f, L, labeling
   !print *, strN, sizeN, nAt, pgOps, diag, a,b,c,d,e,f,L,labeling
   if (status/=0) then; print *,"exiting"; exit; endif
   
   ! Define the full HNF matrix
   HNF = 0
   HNF(1,1) = a; HNF(2,1) = b; HNF(2,2) = c
   HNF(3,1) = d; HNF(3,2) = e; HNF(3,3) = f
   ! Compute the inverse HNF
   call matrix_inverse(real(HNF,dp),Sinv)
   ! Compute the SL inverse
   call matrix_inverse(matmul(p,HNF),sLVinv)

   ! Find the coordinates of the basis atoms
   if(allocated(aBas)) deallocate(aBas)
   allocate(aBas(3,nAt*nD))
   ic = 0
               do i = 1,nD    

   do z1 = 0, a-1
      do z2 = (b*z1)/a, c+(b*z1)/a - 1
         do z3 = z1*(d-(e*b)/c)/a+(e*z2)/c, f+z1*(d-(e*b)/c)/a+(e*z2)/c - 1
               ic = ic + 1; if (ic > nAt*nD) stop "Problem in basis atoms..."
               aBas(:,ic) = matmul(Sinv,(/z1,z2,z3/))+matmul(sLVinv,dvec(:,i))
            enddo
         enddo
      enddo
   enddo
   !print *, ic
   if (ic /= nAt*nD) stop "Not enough basis atoms..."
      write(12,'(a10," str #:",i9)') trim(title), strN
   write(12,'("1.00")')
   ! Make "nice" superlattice vectors (maximally orthogonal, Minkowski reduced)
   sLVorig = matmul(p,HNF)
   call reduce_to_shortest_basis(sLVorig,sLV,eps)
   call matrix_inverse(sLV,sLVinv)
   do i = 1,3
      write(12,'(3f12.8)') sLV(:,i)
   enddo
   do i = 0, k-1
      ic = 0
      do ilab = 1, len(trim(labeling))
         if (labeling(ilab:ilab)==achar(i+48)) ic = ic + 1
      enddo
      write(12,'(i3,1x)',advance='no') ic
   enddo
   write(12,*) ! Start next line
   write(12,'("D")')
   do ilab = 0,k-1
      do i = 1, nAt*nD
         if (labeling(i:i)==achar(ilab+48)) then
            v = matmul(sLVinv,matmul(sLVorig,aBas(:,i)))
            do while(any(v >= 1.0_dp - eps) .or. any(v < 0.0_dp - eps)) 
               v = merge(v, v - 1.0_dp, v <  1.0_dp - eps) 
               v = merge(v, v + 1.0_dp, v >= 0.0_dp - eps)
            enddo
            write(12,'(3f12.8)') v 
         endif
      enddo
   enddo
   write(12,'("0.00 0.00 #Energy and energy factor ",/,"#")')
enddo
END PROGRAM makeStrIn
