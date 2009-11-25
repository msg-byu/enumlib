! This program reads in a file in struct_enum.out format and writes out
! a VASP-like structure file. The program reads in (1) name of the input file,
! and (2) the structure number
PROGRAM makeStr
use num_types
use vector_matrix_utilities
use numerical_utilities
use enumeration_utilities ! This maps structures from real space into the group
implicit none
character(80) fname, title, labeling, strname, strNstring
integer ioerr, iline, ic, i, ilab, pgOps, nD, hnfN, iAt
integer k, strN, istrN, strNi, strNf, sizeN, n, diag(3), a,b,c,d,e,f, HNF(3,3), L(3,3)
integer, pointer :: gIndx(:)
real(dp) :: p(3,3), sLV(3,3), eps, v(3), Ainv(3,3), sLVinv(3,3)
!real(dp) :: map2G(3,3) ! Transformation that takes a real-space lattice point into the group
real(dp), pointer :: aBas(:,:), dvec(:,:) ! pointer so it can be passed in to a routine
character(1) bulksurf
character(40) dummy 

! Testing
integer, pointer :: spin(:) ! Occupation variable of the positions
real(dp), pointer :: x(:) ! Concentration of each component
! Testing

open(13,file="readcheck_makestr.out")
write(13,'(5(/),"This file echoes the input for the makestr.x program and")')
write(13,'("prints some of the computed quantities as well. It is intended to")')
write(13,'("be useful for debugging and fixing bad input")')
write(13,'("***************************************************")')

if(iargc()/=2 .and. iargc()/=3) stop "makestr.x requires two (optional three)  arguments: &
               & the filename to read from and the number of the structure (to start and to stop at)"
call getarg(1,dummy)
read(dummy,*) fname
call getarg(2,dummy)
read(dummy,*) strNi
if (iargc()==3) then
  call getarg(3,dummy)
  read(dummy,*) strNf
else
  strNf=strNi
endif


!call read_nth_line_from_enumlist()

open(11,file=fname,status='old',iostat=ioerr)
if(ioerr/=0)then; write(*,'("Input file doesn''t exist:",a80)') trim(fname);endif
! Read in the title from the struct_enum.out file
read(11,'(a80)') title; title = trim(title)
write(13,'("Title: ",a80)') title

! Read in surf/bulk mode marker
read(11,'(a1)') bulksurf
write(13,'("bulk/surf: ",a1)') bulksurf

! Read in the parent lattice vectors
do i = 1,3; read(11,*) p(:,i); enddo
call matrix_inverse(p,Ainv)
write(13,'("Parent lattice:",/,3(3f7.3,1x,/))') (p(:,i),i=1,3) 
write(13,'("Parent lattice inverse matrix:",/,3(3f7.3,1x,/))') (Ainv(:,i),i=1,3) 

! Read in the number of d-vectors, then read each one in
read(11,*) nD
allocate(dvec(3,nD))
do i = 1,nD; read(11,*) dvec(:,i); enddo
write(13,'("d-vectors:",/,40(3f8.4,1x,/))') (dvec(:,i),i=1,nD) 


! Read in the number of labels, i.e., binary, ternary, etc.
read(11,'(i2)') k
write(13,'("Number of labels (k): ",i2)') k 
read(11,*)
read(11,*) eps
write(13,'("Finite precision parameter (epsilon): ",g17.10)') eps 

do 
   read(11,*) dummy
   if(dummy=="start") exit
enddo
write(13,'("Found ""start"" label. Now skipping lines to struct #: ",i10)') strNi 

! Skip to the specified structure number
do iline = 1, strNi-1
   read(11,*)! dummy; print *,dummy
enddo
write(13,'("Skipped to the ",i10,"nth line in ",a80)') strNi,fname 

!print *,"<<<<< makestr.x should be updated to call map_enumStr_to_real_space >>>>>"
! Read in the info for the given structure
do istrN=strNi,strNf
   read(11,*) strN, hnfN, sizeN, n, pgOps, diag, a,b,c,d,e,f, L, labeling
   write(13,'("strN, hnfN, sizeN, n, pgOps, diag, a,b,c,d,e,f, L, labeling")')  
   write(13,'(i11,1x,i7,1x,i8,1x,i2,1x,i2,1x,3(i3,1x),1x,6(i3,1x),1x,9(i3,1x),1x,a40)') &
        strN, hnfN, sizeN, n, pgOps, diag, a,b,c,d,e,f, L, labeling
   
   L = transpose(L) ! Listed with columns as the fast index in the struct_enum.out but reads
                    ! in with rows as the fast index. This was the source of a hard-to-find bug

   ! Define the full HNF matrix
   HNF = 0
   HNF(1,1) = a; HNF(2,1) = b; HNF(2,2) = c
   HNF(3,1) = d; HNF(3,2) = e; HNF(3,3) = f
   

   call map_enumStr_to_real_space(k,n,HNF,labeling,p,dvec,eps,sLV,aBas,spin,gIndx,x,L,diag)

!   gIndx = -1
!   
!   write(13,'("HNF matrix:",/,3(3(i2,1x),/))') (HNF(i,:),i=1,3) 
!   ! Compute the inverse HNF 
!   call matrix_inverse(real(HNF,dp),Sinv)
!   write(13,'("Invers HNF matrix:",/,3(3(f7.3,1x),/))') (Sinv(i,:),i=1,3) 
!   
!   ! Compute the superlattice vectors
!   sLV = matmul(p,HNF) 
!   call matrix_inverse(sLV,sLVinv)
!   write(13,'("Superlattice inv vectors:",/,3(3(f7.3,1x),/))') (sLVinv(i,:),i=1,3) 
!   
!   ! Compute the mapping for x -> G
!   !print *,"shape of L*Ainv",shape(matmul(L,Ainv))
!   !print *,"shape of diag",shape(diag)
!   
!   map2G = matmul(L,Ainv)
!   write(13,'("x->G mapping:",/,3(3(f7.3,1x),/))') (map2G(i,:),i=1,3) 
!   
!   !do i=1,3; map2G(:,i)=modulo(map2G(:,i),real(diag,dp));enddo
!   !write(13,'("x->G mapping:",/,3(3(f7.3,1x),/))') (map2G(i,:),i=1,3) 
!   
!   write(13,'("Coordinates of the atoms:")')  
!   ! Find the coordinates of the basis atoms
!   allocate(aBas(3,nAt*nD))
!   ic = 0
!   do i = 1,nD
!   do z1 = 0, a-1
!      do z2 = (b*z1)/a, c+(b*z1)/a - 1
!         do z3 = z1*(d-(e*b)/c)/a+(e*z2)/c, f+z1*(d-(e*b)/c)/a+(e*z2)/c - 1
!               ic = ic + 1; if (ic > nAt*nD) stop "Problem in basis atoms..."
!               ! Matmul by Sinv puts the position of the parent unit cell into Cartesian coordinates
!               ! and adding the d-vector on just gives the shift inside the current parent cell
!               aBas(:,ic) = matmul(Sinv,(/z1,z2,z3/))+matmul(sLVinv,dvec(:,i))
!               ! This matmul takes the Cartesian coordinates of the vector and converts them to direct
!               ! (lattice) coordinates of the supercell
!               aBas(:,ic) = matmul(sLV,aBas(:,ic))
!               write(13,'("at #: ",i3,4x,"position: ",3(f7.3,1x))') ic,aBas(:,ic) 
!   
!               ! Now take aBas and map it into the group.
!               ! (LA^-1)^-1 takes a real space vector and maps it into the group.
!               ! But subtract the d-set member since we only care about which parent unit cell this
!               ! atom is in
!               greal = matmul(map2G,matmul(sLV,matmul(Sinv,(/z1,z2,z3/))) )
!               g = nint(greal)
!               write(13,'("group member (real):",3(f8.3,1x))') greal 
!               write(13,'("group member (nint):",3(i3,1x))') g
!               if(.not. equal(greal,g,eps)) stop "map2G didn't work"
!               g = modulo(g,diag)
!               write(13,'("group member modded:",3(i3,1x))') g
!               ! Map the group element components (and d-vector index) to a single number
!               ! This indexes the group member in the labeling
!               gIndx(ic) = (i-1)*nAt*nD+g(1)*diag(2)*diag(3) + g(2)*diag(3) + g(3) + 1
!               write(13,'("Ordinal # in the group:",i3)') gIndx(ic)
!               !write(*,'(3f12.8,5x,3f12.8)') matmul(Sinv,(/z1,z2,z3/)),matmul(sLVinv,dvec(:,i))
!            enddo
!         enddo
!      enddo
!   enddo
!   if (ic /= nAt*nD) stop "Not enough basis atoms..."
!   write(13,'("Atom labels and positions:",/,3(3(f7.3,1x),/))')  
!   do iAt = 1, nAt*nD
!      write(13,'("gIndx: ",i2," Atom #: ",i2,1x,"label: ",a1,1x,3(f7.3,1x))') &
!           gIndx(iAt), iAt, labeling(gIndx(iAt):gIndx(iAt)), aBas(:,iAt)
!   enddo
   

   write(strname,'("vasp.",i4.4)') strN
   write(strNstring,*) strN
   open(12,file=strname)
   write(12,'(a80)') trim(adjustl(title)) // " str #: " // adjustl(strNstring)
   write(12,'("scale factor")')
   sLV=matmul(p,HNF)
   call matrix_inverse(sLV,sLVinv)
   do i = 1,3
      write(12,'(3f12.8)') sLV(:,i)
   enddo
   write(13,'("New inverse after reduction",/,3(3(f7.3,1x),/))') (sLVinv(i,:),i=1,3) 
   
   ! This part counts the number of atoms of each type and lists the numbers on onlien before the atomic basis vectors
   do i = 0, k-1
      ic = 0
      do iAt = 1, n*nD
         if (labeling(gIndx(iAt):gIndx(iAt))==char(i+48)) then
            ic = ic + 1
            write(13,'("Atom #: ",i2," given label #",i2," label: ",a1)') iAt, ic,labeling(gIndx(iAt):gIndx(iAt))
         endif
      enddo
      write(12,'(i3,1x)',advance='no') ic
      write(13,'(i3," atoms of label type ",i2," found")') ic, i
   enddo
   write(13,'("Finished counting the atoms of each type")')  
   
   write(12,*) ! Start next line
   write(12,'("D")')
   
   ! This part lists the atomic basis vectors that we found in the triple z1, z2, z3 loops above.
   ! For vasp, UNCLE it needs to list the vectors in blocks of that have the same label.
   call cartesian2direct(sLV,aBas,eps)
   do ilab = 0,k-1
      do iAt = 1, n*nD
         if (labeling(gIndx(iAt):gIndx(iAt))==char(ilab+48)) then 
!!            v = matmul(sLVinv,aBas(:,iAt)) ! Put positions into "direct" coordinates
!!            ! This keeps the atomic coordinates inside the first unit cell---
!!            ! not necessary but aesthetically pleasing.
!!            do while(any(v >= 1.0_dp - eps) .or. any(v < 0.0_dp - eps)) 
!!               v = merge(v, v - 1.0_dp, v <  1.0_dp - eps) 
!!               v = merge(v, v + 1.0_dp, v >= 0.0_dp - eps)
!!            enddo
            write(12,'(3f12.8)') aBas(:,iAt)
            write(13,'("At. #: ",i2," position:",3(f7.3,1x),"<",a1,">")') iAt, v, labeling(gIndx(iAt):gIndx(iAt))
            !write(*,'(3f12.8)') aBas(:,i)
         endif
      enddo
   enddo
   close(12);
      
!   do iAt = 1,nAt*nD
!      write(14,'("At#: ",i2," pos: ",3(f7.3,1x),"gIndx: ",i1," label",a1)') iAt, aBas(:,iAt),&
!           & gIndx(iAt), labeling(gIndx(iAt):gIndx(iAt))
!   enddo
!
   deallocate(spin,x,gIndx,aBas)
enddo ! structure loop
close(11)
close(13)
END PROGRAM makeStr
