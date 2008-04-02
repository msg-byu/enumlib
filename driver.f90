PROGRAM driver
use num_types 
use derivative_structure_generator
use labeling_related
use rational_mathematics
use crystal_types
use io_utilities
implicit none

integer ivol, ivolTot, iuq, ihnf, i, ilab ! Loop counters
integer Nq, nMin, nMax, runTot  ! Numbers of various things
integer diag(3), ld(6) ! diagonal elements of SNF, lower diagonal---elements of HNF
integer k ! Number of colors/label types (i.e., binary, ternary, etc.)
integer LatDim ! 2D or 3D parent lattice?
integer, pointer :: labelings(:,:) =>null(), SNF_labels(:) =>null(), tlab(:,:)=>null()
integer, pointer :: HNF(:,:,:) => null(), reducedHNF(:,:,:) => null(), G(:,:)
integer, pointer :: trgroup(:,:)=>null()
real(dp), pointer :: uqlatts(:,:,:) => null()
integer, pointer, dimension(:,:,:) :: SNF => null(), L => null(), B => null()
integer, pointer, dimension(:,:,:) :: uqSNF => null()
integer :: ctot, csize ! counters for total number of structures and structures of each size
real(dp) tend, tstart, eps
real(dp) :: parLV(3,3)
type(opList), pointer :: fixOp(:)
character, pointer :: table(:)
character(80) title

call read_input(title,LatDim,parLV,k,nMin,nMax,eps) ! Read in parent lattice vectors, etc.
title = trim(title)
open(14,file="struct_enum.out")
write(14,'(a80)') title
do i = 1,3
   write(14,'(3(e14.8,1x),3x,"# a",i1," parent lattice vector")') parLV(:,i),i
enddo
write(14,'(i2,"-nary case")') k
write(14,'(8x,"#tot",5x,"#size",1x,"nAt",2x,"pg",4x,"SNF",13x,"HNF",17x,"Left transform",17x,"labeling")')
 
write(*,'("Volume",7x,"CPU",5x,"#HNFs",3x,"#SNFs",&
          &4x,"#reduced",4x,"% dups",6x,"volTot",6x,"RunTot")')
!allocate(tempStr1(1))
runTot = 0
ctot = 0
do ivol = max(k,nMin),nMax
   csize = 0
   call cpu_time(tstart)
   if (LatDim==3) then !<<< 2D ? or 3D? >>
      call get_all_HNFs(ivol,HNF)  ! 3D
   else
      call get_all_2D_HNFs(ivol,HNF) ! 2D
   endif
   call remove_duplicate_lattices(HNF,LatDim,parLV,reducedHNF,fixOp,uqlatts,eps)
   call get_SNF(reducedHNF,L,SNF,B,uqSNF,SNF_labels)
   ! We have a list containing each unique derivative superlattice and
   ! its corresponding Smith Normal Form. Now make the labelings.
   Nq = size(uqSNF,3)  ! Number of unique SNFs
   ivolTot = 0
   do iuq = 1, Nq ! Loop over all of the unique SNFs
      diag = (/uqSNF(1,1,iuq),uqSNF(2,2,iuq),uqSNF(3,3,iuq)/)
      call make_member_list(diag,G)  ! Need the image group G for removing lab-rot dups
      call generate_labelings(k,diag,labelings,table,trgroup) ! Removes trans-dups,non-prims,label-exchange dups
      do ihnf = 1, size(SNF_labels) ! Store each of the HNF and left transformation matrices
         if (SNF_labels(ihnf)/=iuq) cycle ! Skip structures that don't have the current SNF
         ! Need to do this step for each HNF, not each SNF
         allocate(tlab(size(labelings,1),size(labelings,2))); tlab = labelings
         ! tlab is a temporary copy of the labelings (only the unique ones)
         ! table is a list of markers (D,F,E,etc.) for duplicate labelings
         ! G is a collection of 3-vectors representing members of the quotient group
         ! trgroup is a list of all the permutations that correspond to translations
         call remove_label_rotation_dups(L(:,:,iHNF),parLV,fixOp(iHNF)%rot,G,&
                                        tlab,table,trgroup,k,diag,eps)
         ivolTot = ivolTot + size(tlab,1)
         ld = (/reducedHNF(1,1,ihnf),reducedHNF(2,1,ihnf),reducedHNF(2,2,ihnf),&
                reducedHNF(3,1,ihnf),reducedHNF(3,2,ihnf),reducedHNF(3,3,ihnf)/)
         do ilab = 1,size(tlab,1) ! write out the labelings to a file
            ctot = ctot + 1
            csize = csize + 1
            write(14,'(i11,1x,i9,1x,i3,2x,i3,2x,3(i2,1x),2x,6(i2,1x),2x,9(i4),2x,40i1)') &
              ctot, csize,ivol,size(fixOp(ihnf)%rot,3),diag,ld,&
              transpose(L(:,:,iHNF)),tlab(ilab,:)
         enddo
      enddo
   enddo
   call cpu_time(tend)
   runTot = runTot + ivolTot
   write(*,'(i4,1x,f12.2,1x,i8,3x,i3,3x,i7,7x,f7.4,i12,i12)')&
        ivol,tend-tstart,size(HNF,3),&
        size(uqSNF,3),size(reducedHNF,3),1-size(reducedHNF,3)/real(size(HNF,3)),&
                ivolTot, runTot
enddo
close(14)

END PROGRAM driver
