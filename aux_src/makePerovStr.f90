! This program reads in a file in struct_enum.out format and writes out
! a VASP-like structure file. The program reads in (1) name of the input file,
! and (2) the structure number
PROGRAM makeStr
use num_types
use vector_matrix_utilities
use cwrapper
implicit none
character(80) fname, title, strname, strNstring
character(100) :: labeling
integer ioerr, iline, z1, z2, z3, ic, i, ilab, pgOps, nD, hnfN, iFace
integer k, strN, sizeN, nAt, diag(3), a,b,c,d,e,f, HNF(3,3), L(3,3)
real(dp) :: p(3,3), sLV(3,3), Sinv(3,3), sLVorig(3,3), eps, v(3), Binv(3,3), sLVinv(3,3), face(3,3)
real(dp), allocatable :: aBas(:,:), dvec(:,:)
character(1) bulksurf
character(40) dummy 
character(2) label(4)

face(:,1) = (/0.5_dp,0.5_dp,0._dp /)
face(:,2) = (/0.5_dp,0._dp ,0.5_dp/)
face(:,3) = (/0._dp ,0.5_dp,0.5_dp/)
label(1) = "A"
label(2) = "A'"
label(3) = "B"
label(4) = "Ox"
if(iargc()/=2) stop "makeperovstr.x requires two arguments: &
               & the filename to read from and the number of the structure"
call getarg(1,dummy)
read(dummy,*) fname
call getarg(2,dummy)
read(dummy,*) strN

!call read_nth_line_from_enumlist()

open(11,file=fname,status='old',iostat=ioerr)
if(ioerr/=0)then; write(*,'("Input file doesn''t exist:",a80)') trim(fname);endif
! Read in the title from the struct_enum.out file
read(11,'(a80)') title; title = trim(title)

! Read in surf/bulk mode marker
read(11,'(a1)') bulksurf

! Read in the parent lattice vectors
do i = 1,3; read(11,*) p(:,i); enddo
call matrix_inverse(p,Binv)

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
read(11,*) strN, hnfN, sizeN, nAt, pgOps, diag, a,b,c,d,e,f, L, labeling
close(11)
!print *, strN, sizeN, nAt, pgOps, diag, a,b,c,d,e,f,L,labeling

! Define the full HNF matrix
HNF = 0
HNF(1,1) = a; HNF(2,1) = b; HNF(2,2) = c
HNF(3,1) = d; HNF(3,2) = e; HNF(3,3) = f
! Compute the inverse HNF 
call matrix_inverse(real(HNF,dp),Sinv)
! Compute the superlattice vectors 
call matrix_inverse(matmul(p,HNF),sLVinv)

! Find the coordinates of the basis atoms
allocate(aBas(3,nAt*nD))
ic = 0
do i = 1,nD
do z1 = 0, a-1
   do z2 = (b*z1)/a, c+(b*z1)/a - 1
      do z3 = z1*(d-(e*b)/c)/a+(e*z2)/c, f+z1*(d-(e*b)/c)/a+(e*z2)/c - 1
            ic = ic + 1; if (ic > nAt*nD) stop "Problem in basis atoms..."
            aBas(:,ic) = matmul(Sinv,(/z1,z2,z3/))+matmul(sLVinv,dvec(:,i))
            !write(*,'(3f12.8,5x,3f12.8)') matmul(Sinv,(/z1,z2,z3/)),matmul(sLVinv,dvec(:,i))
         enddo
      enddo
   enddo
enddo
if (ic /= nAt*nD) stop "Not enough basis atoms..."

write(strname,'("vasp.",i4.4)') strN
write(strNstring,*) strN
open(12,file=strname)
write(12,'(a80)') trim(adjustl(title)) // " str #: " // adjustl(strNstring)
write(12,'("scale factor")')
! Make "nice" superlattice vectors (maximally orthogonal, Minkowski reduced)
sLVorig = matmul(p,HNF)
call aflow_reduce_to_shortest_basis(sLVorig,sLV,eps)
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
! Write out the "unalloyed non oxygen atoms" and the oxygens
write(12,'(3(i3,1x))',advance='no') nAt*nD, nAt*nD*3
write(12,*) ! Start next line
write(12,'("D")')

do ilab = 0,k-1 ! Write out the "alloy" atoms 
   do i = 1, nAt*nD
      if (labeling(i:i)==achar(ilab+48)) then 
         v = matmul(sLVinv,matmul(sLVorig,aBas(:,i)))
         do while(any(v >= 1.0_dp - eps) .or. any(v < 0.0_dp - eps)) 
            v = merge(v, v - 1.0_dp, v <  1.0_dp - eps) 
            v = merge(v, v + 1.0_dp, v >= 0.0_dp - eps)
         enddo
         write(12,'(3f12.8,1x,a2)') v,label(ilab+1)
         !write(*,'(3f12.8)') aBas(:,i)
      endif
   enddo
enddo

do i = 1, nAt*nD ! Write out the B atoms
   v = matmul(sLVinv,matmul(sLVorig,aBas(:,i))+(/.5_dp,.5_dp,.5_dp/))
   do while(any(v >= 1.0_dp - eps) .or. any(v < 0.0_dp - eps)) 
      v = merge(v, v - 1.0_dp, v <  1.0_dp - eps) 
      v = merge(v, v + 1.0_dp, v >= 0.0_dp - eps)
   enddo
   write(12,'(3f12.8,1x,a2)') v, label(3)
enddo

do i = 1, nAt*nD ! Write out the oxygen atoms
   do iFace = 1,3
      v = matmul(sLVinv,matmul(sLVorig,aBas(:,i))+face(:,iFace))
      do while(any(v >= 1.0_dp - eps) .or. any(v < 0.0_dp - eps)) 
         v = merge(v, v - 1.0_dp, v <  1.0_dp - eps) 
         v = merge(v, v + 1.0_dp, v >= 0.0_dp - eps)
      enddo
      write(12,'(3f12.8,1x,a2)') v, label(4)
   enddo
enddo


END PROGRAM makeStr
