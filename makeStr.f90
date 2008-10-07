! This program reads in a file in struct_enum.out format and writes out
! a VASP-like structure file. The program reads in (1) name of the input file,
! and (2) the structure number
PROGRAM makeStr
use num_types
use vector_matrix_utilities
implicit none
character(80) fname, title, labeling, strname
integer ioerr, iline, z1, z2, z3, ic, i, ilab, pgOps, nD
integer k, strN, sizeN, nAt, diag(3), a,b,c,d,e,f, HNF(3,3), L(3,3)
real(dp) :: p(3,3), sLV(3,3), Sinv(3,3), sLVorig(3,3), eps
real(dp), allocatable :: aBas(:,:), dvec(:,:)
character(1) bulksurf
!character(40) dummy 
read(*,*) fname, strN
!fname = "out.fcc.binary.17"; strN=14
open(11,file=fname,status='old',iostat=ioerr)
if(ioerr/=0)then; write(*,'("Input file doesn''t exist:",a80)') trim(fname);endif
! Read in the title from the struct_enum.out file
read(11,'(a80)') title; title = trim(title)

! Reand in surf/bulk mode marker
read(11,'(a1)') bulksurf

! Read in the parent lattice vectors
do i = 1,3; read(11,*) p(:,i); enddo

! Read in the number of d-vectors, then read each one in
read(11,*) nD
allocate(dvec(3,nD))
do i = 1,nD; read(11,*) dvec(:,i); enddo

! Read in the number of labels, i.e., binary, ternary, etc.
read(11,'(i2)') k
read(11,*)
read(11,*) eps

! Skip to the specified structure number
do iline = 1, strN+1
   read(11,*)! dummy; print *,dummy
enddo
! Read in the info for the given structure
read(11,*) strN, sizeN, nAt, pgOps, diag, a,b,c,d,e,f, L, labeling
close(11)
!print *, strN, sizeN, nAt, pgOps, diag, a,b,c,d,e,f,L,labeling

! Define the full HNF matrix
HNF = 0
HNF(1,1) = a; HNF(2,1) = b; HNF(2,2) = c
HNF(3,1) = d; HNF(3,2) = e; HNF(3,3) = f
! Compute the superlattice vectors 
sLVorig = matmul(p,HNF)
call matrix_inverse(real(HNF,dp),Sinv)
! Make "nice" superlattice vectors (maximally orthogonal, Minkowski reduced)
call reduce_to_shortest_basis(sLVorig,sLV,eps)

! Find the coordinates of the basis atoms
allocate(aBas(3,nAt*nD))
ic = 0
do i = 1,nD
do z1 = 0, a-1
   do z2 = (b*z1)/a, c+(b*z1)/a - 1
      do z3 = z1*(d-(e*b)/c)/a+(e*z2)/c, f+z1*(d-(e*b)/c)/a+(e*z2)/c - 1
         ic = ic + 1; if (ic > nAt*nD) stop "Problem in basis atoms..."
         aBas(:,ic) = matmul(Sinv,(/z1,z2,z3/))+matmul(Sinv,dvec(:,i))
      enddo
   enddo
enddo
enddo
if (ic /= nAt*nD) stop "Not enough basis atoms..."

write(strname,'("vasp.",i4.4)') strN
open(12,file=strname)
write(12,'(a10," str #:",i9)') trim(title), strN
write(12,'("scale factor")')
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
      if (labeling(i:i)==achar(ilab+48)) write(12,'(3f12.8)') aBas(:,i)
   enddo
enddo

END PROGRAM makeStr
