! Generate all derivative structures of a parent structure
! Gus L. W. Hart BYU July 2007
! Things still to do:
! * Currently only applies to primitive parent structures*

MODULE derivative_structure_generator
use num_types
use enumeration_types
use labeling_related
use combinatorics
use symmetry_module
use compare_structures
use rational_mathematics
use numerical_utilities
use vector_matrix_utilities, only: determinant, matrix_inverse

implicit none
private
public get_all_HNFs, remove_duplicate_lattices, get_SNF, get_all_2D_HNFs, generate_derivative_structures
CONTAINS

!***************************************************************************************************
! Finds all the possible diagonals of the HNF matrices of a given size
SUBROUTINE get_HNF_diagonals(detS, diagonals)
integer, intent(in) :: detS ! Cell size, i.e., determinant of S matrix
integer, pointer :: diagonals(:,:) ! All possible diagonals

integer i, j, id, quotient, status
integer :: tempDiag(3,detS*3)

id = 0 ! Number of diagonals found
do i = 1,detS ! Loop over possible first factors
   if (.not. mod(detS,i)==0) cycle
   quotient = detS/i
   do j = 1,quotient  ! Loop over possible second/third factors
      if (.not. mod(quotient,j)==0) cycle
      id = id + 1
      tempDiag(:,id) = (/i, j, quotient/j /) ! Construct the factor triplet
   enddo
enddo
allocate(diagonals(3,id),STAT=status)
if(status/=0) stop "Allocation failed in get_HNF_diagonals, module deriv..."
diagonals = tempDiag(:,1:id)
END SUBROUTINE get_HNF_diagonals

!*******************************************************************************
! This subroutine generates all the unique HNF matrices of a given determinant.
! (See Santoro and Mighell 1972 Acta. Cryst.)
SUBROUTINE get_all_HNFs(volume,hnf)
integer, intent(in) :: volume
integer, pointer:: hnf(:,:,:)

integer, pointer    :: d(:,:) => null()
integer             :: i, j, k, l    ! Loop counters
integer             :: N, Nhnf, ihnf ! # of triplets, # of HNF matrices, HNF counter
integer status
call get_HNF_diagonals(volume,d)
N = size(d,2)

! Count the total number of HNF matrices for given determinant (volume)
Nhnf = 0
do i = 1,N
   Nhnf = Nhnf + d(2,i)*d(3,i)**2
enddo

allocate(hnf(3,3,Nhnf),STAT=status)
if(status/=0) stop "Failed to allocate memory in get_all_HNFs"
ihnf = 0
do i = 1,N ! Loop over the permutations of the diagonal elements of the HFNs
   do j = 0,d(2,i)-1  ! Look over possible values of row 2, element 1
   do k = 0,d(3,i)-1  ! Ditto for row 3, element 1
   do l = 0,d(3,i)-1  ! Ditto for row 3, element 2
      ihnf = ihnf+1 ! Count the HNFs and construct the next one
      hnf(:,:,ihnf) = reshape((/ d(1,i),      0,     0,        &   
                                      j, d(2,i),     0,        &   
                                      k,      l, d(3,i)  /), (/3,3/))
      hnf(:,:,ihnf) = transpose(hnf(:,:,ihnf))
   enddo;enddo;enddo  ! End loops over values for off-diagonal elements
enddo ! End loop over all unique triplets of target determinant (volume)
if (ihnf /= Nhnf) stop "HNF: not all the matrices were generated...(bug!)"
END SUBROUTINE get_all_HNFs

!*******************************************************************************
! Find all the SNFs of a group of HNFs. A and B are the transformation matrices.
! F is a list of the *unique* SNFs. SNF_label indicates which of the unique SNFs
! corresponds to each HNF.
SUBROUTINE get_SNF(HNF,A,SNF,B,F,SNF_label)
integer, pointer :: HNF(:,:,:), SNF_label(:)
integer, pointer, dimension(:,:,:) :: A, SNF, B, F

integer ihnf, nHNF, nfound, ifound, status
logical duplicate
integer, allocatable, dimension(:,:,:) :: tF

nHNF = size(HNF,3)
nfound = 0
allocate(A(3,3,nHNF),B(3,3,nHNF),SNF(3,3,nHNF),tF(3,3,nHNF),SNF_label(nHNF),STAT=status)
if(status/=0) stop "Failed to allocate memory in get_SNF"

do ihnf = 1, nHNF ! Loop over each HNF in the list
   call SmithNormalForm(HNF(:,:,ihnf),A(:,:,ihnf),SNF(:,:,ihnf),B(:,:,ihnf))
   duplicate = .false.
   do ifound = 1,nfound ! Check to see if this SNF is unique or a duplicate.
      if (all(SNF(:,:,ihnf)==tF(:,:,ifound))) then
         duplicate = .true.; SNF_label(ihnf) = ifound; exit; endif
   enddo
   if (.not. duplicate) then ! store this SNF in the unique list
      nfound = nfound + 1
      tF(:,:,nfound) = SNF(:,:,ihnf)
      SNF_label(ihnf) = nfound
   endif
enddo
if (associated(F)) deallocate(F)
allocate(F(3,3,nfound),STAT=status)
if(status/=0) stop "Failed to allocate memory in get_SNF for array F"
F = tF(:,:,1:nfound)

ENDSUBROUTINE get_SNF

!*******************************************************************************
! Takes a bunch of HNF matrices and a parent lattice and removes those that are
! rotationally equivalent (under the rotations of the parent lattice). Also
! returns all the unique derivative _lattices_ for this parent. Returns also 
! a list of rotations that fix each superlattice. 
SUBROUTINE remove_duplicate_lattices(hnf,LatDim,parent_lattice,uq_hnf,fixing_op,latts,eps)
integer, pointer :: hnf(:,:,:) ! HNF matrices (input)
integer :: LatDim ! Is the parent lattice 2D or 3D?
real(dp), intent(in) :: parent_lattice(3,3) ! parent lattice (input)
integer, pointer :: uq_hnf(:,:,:) ! list of symmetrically distinct HNFs (output)
type(opList), pointer :: fixing_op(:) ! List of operations that fix each HNF
real(dp),intent(in):: eps ! finite precision (input)

real(dp), pointer:: sgrots(:,:,:)
real(dp), dimension(3,3) :: test_latticei, test_latticej, thisRot, rotLat, origLat
integer i, Nhnf, iuq, irot, j, nRot, Nq, ic, status
integer, allocatable :: temp_hnf(:,:,:)
real(dp), pointer :: latts(:,:,:)
logical duplicate
type(opList), allocatable :: tmpOp(:)
real(dp), allocatable :: tSGrots(:,:,:)

!do i=1,3
! write(*,'(3f10.6)') parent_lattice(:,i)
!enddo
call get_lattice_pointGroup(parent_lattice,sgrots,eps)
nRot = size(sgrots,3)
!print *,"<<< parent Ops found:",nRot

Nhnf = size(hnf,3)
allocate(temp_hnf(3,3,Nhnf),STAT=status)
if(status/=0) stop "Failed to allocate memory in remove_duplicate_lattices"
temp_hnf = hnf

! For the 2D case, eliminate the "3D" operations.
if (LatDim==2) then
   allocate(tSGrots(3,3,nRot),STAT=status)
   if(status/=0) stop "Failed to allocate memory in remove_duplicate_lattices: tSGrots"
   irot = 0
   do i = 1, nRot
      if (equal(sgrots(2:3,1,i),0._dp,eps) .and. equal(sgrots(1,2:3,i),0._dp,eps) .and.&
           & equal(sgrots(1,1,i),1._dp,eps)) then ! this operation is "2D"
         irot = irot + 1
         tSGrots(:,:,irot) = sgrots(:,:,i)
         !do j = 1,3  ! Take this out after debugging...
         !   print *, irot, sgrots(j,:,i)
         !   !write(10,'(3f8.3)') sgrots(j,:,i)
         !enddo
         !print *
      endif
   enddo
   nRot = irot
   deallocate(sgrots)
   allocate(sgrots(3,3,nRot),STAT=status)
if(status/=0) stop "Allocation of sgrots failed in remove_duplicate_lattices, module deriv..."
   sgrots=tSGrots(:,:,1:nRot)
endif

! For each HNF in the list, see if it is a derivative lattice of a preceding
! HNF in the list. If so, don't include it in the list of unique ones.
iuq = 1
do i = 2,Nhnf  ! Loop over each matrix in the original list
   duplicate = .false.
   do j = 1,iuq! Loop over the known unique matrices (in the updated list)
      do irot = 1, nRot! The duplicates will always be rotated from 
                       ! the original (otherwise necessarily unique)
         test_latticei = matmul(sgrots(:,:,irot),matmul(parent_lattice,hnf(:,:,i)))
         test_latticej = matmul(parent_lattice,temp_hnf(:,:,j))
         if (is_equiv_lattice(test_latticei,test_latticej,eps)) then
            duplicate = .true.;exit;endif
      enddo
   enddo
   if (.not. duplicate) then ! No duplicate found---add this to the unique list
      iuq = iuq+1
      temp_hnf(:,:,iuq) = hnf(:,:,i);endif
enddo
allocate(uq_hnf(3,3,iuq),latts(3,3,iuq),STAT=status)
if(status/=0) stop "Failed to allocate uq_hnf/latts in remove_duplicate_lattices: uq_hnf"
uq_hnf = temp_hnf(:,:,:iuq)
Nq = iuq
forall (i=1:Nq); latts(:,:,i) = matmul(parent_lattice,uq_hnf(:,:,i));end forall
! Now loop through the list of unique superlattices and record which of the parent lattice
! symmetry operations leave the superlattice unchanged. These operations, while they leave
! the superlattice unchanged can permute the labels inside. Thus, this is another source of 
! duplicates. 
allocate(tmpOp(Nq),fixing_op(Nq),STAT=status)
if(status/=0) stop "Allocation of tmpOp or fixing_op failed in remove_duplicate_lattices, module deriv..."
do i=1,Nq; allocate(tmpOp(i)%rot(3,3,nRot),STAT=status)
if(status/=0) stop "Allocation failed in get_HNF_diagonals, module deriv...";enddo
do iuq = 1, Nq  ! Loop over each unique HNF
   ic = 0
   do iRot = 1,nRot  ! Loop over each rotation
      thisRot = sgrots(:,:,iRot) ! Store the rotation
      origLat = matmul(parent_lattice,uq_hnf(:,:,iuq))  ! Compute the superlattice
      rotLat = matmul(thisRot,origLat)          ! Compute the rotated superlattice
      if (is_equiv_lattice(rotLat,origLat,eps)) then ! this operation fixes the lattice and should be recorded
         ic = ic + 1
         tmpOp(iuq)%rot(:,:,ic) = thisRot
      endif
   enddo ! Now we know which rotations fix the lattice and how many there are so store them
   do i=1,ic; allocate(fixing_op(iuq)%rot(3,3,ic),STAT=status)
if(status/=0) stop "Allocation of fixing_op%rot failed, module deriv..."; enddo ! Allocate the storage for them
   fixing_op(iuq)%rot = tmpOp(iuq)%rot(:,:,1:ic) ! Stuff the rotations into the permanent array
enddo

END SUBROUTINE remove_duplicate_lattices

!*******************************************************************************
! This subroutine generates all the unique HNF matrices of a given determinant
! (See Santoro and Mighell 1972 Acta. Cryst.) but for a quasi-2D case
SUBROUTINE get_all_2D_HNFs(volume,hnf)
integer, intent(in) :: volume
integer, pointer:: hnf(:,:,:)

integer, pointer    :: d(:,:) => null()
integer             :: i, j    ! Loop counters
integer             :: N, Nhnf, ihnf ! # of triplets, # of HNF matrices, HNF counter
integer status
call get_HNF_2D_diagonals(volume,d)
N = size(d,2)

! Count the total number of HNF matrices for given determinant (volume)
Nhnf = 0
do i = 1,N
   Nhnf = Nhnf + d(3,i)
enddo

allocate(hnf(3,3,Nhnf),STAT=status)
if(status/=0) stop "Failed to allocate memory in get_all_2D_HNFs"
ihnf = 0
do i = 1,N ! Loop over the permutations of the diagonal elements of the HFNs
   do j = 0,d(3,i)-1  ! Look over possible values of row 2, element 1
      ihnf = ihnf+1 ! Count the HNFs and construct the next one
      hnf(:,:,ihnf) = reshape((/ d(1,i),      0,     0,        &   
                                      0, d(2,i),     0,        &   
                                      0,      j, d(3,i)  /), (/3,3/))
      hnf(:,:,ihnf) = transpose(hnf(:,:,ihnf))
   enddo  ! End loops over values for off-diagonal elements
enddo ! End loop over all unique triplets of target determinant (volume)
if (ihnf /= Nhnf) stop "HNF: not all the matrices were generated...(bug!)"
END SUBROUTINE get_all_2D_HNFs

!***************************************************************************************************
! Finds all the possible diagonals of the HNF matrices of a given size
SUBROUTINE get_HNF_2D_diagonals(detS, diagonals)
integer, intent(in) :: detS ! Cell size, i.e., determinant of S matrix
integer, pointer :: diagonals(:,:) ! All possible diagonals

integer i, id, quotient, status
integer :: tempDiag(3,detS*3)

id = 0 ! Number of diagonals found
do i = 1,detS ! Loop over possible first factors
   if (.not. mod(detS,i)==0) cycle
   quotient = detS/i
   id = id + 1
   tempDiag(:,id) = (/1, i, quotient /) ! Construct the factor triplet
enddo
allocate(diagonals(3,id),STAT=status)
if(status/=0) stop "Allocation of diagonals fails in get_HNF_2D_diagonals, module deriv..."
diagonals = tempDiag(:,1:id)
END SUBROUTINE get_HNF_2D_diagonals

!****************************************************************************************************
SUBROUTINE generate_derivative_structures(title, parLV, k, nMin, nMax, pLatTyp, eps, full)
integer, intent(in) :: k, nMin, nMax 
character(10), intent(in) :: title
real(dp), intent(in) :: parLV(3,3)
character(1), intent(in) :: pLatTyp
logical, intent(in) :: full 

integer, pointer, dimension(:,:,:) :: uqSNF => null()
integer, pointer :: trgroup(:,:)=>null()
real(dp), pointer :: uqlatts(:,:,:) => null()
character, pointer :: table(:)
integer i, ihnf, ilab, iuq, ivol, Nq, iVolTot, runTot, LatDim, status, ilr, nHNF, nRedHNF
integer irg, nrg
integer, pointer :: HNF(:,:,:) => null(), reducedHNF(:,:,:) => null(), G(:,:) => null()
integer, pointer :: labelings(:,:) =>null(), SNF_labels(:) =>null(), tlab(:,:)=>null()
integer, pointer, dimension(:,:,:) :: SNF => null(), L => null(), B => null()
type(opList), pointer :: fixOp(:)
integer :: ctot, csize ! counters for total number of structures and structures of each size
integer diag(3), ld(6) ! diagonal elements of SNF, lower diagonal---elements of HNF
real(dp) eps, tstart, tend, tHNFs, tDupLat, tSNF, tiuq, tML, tGenLab, tRotDup
integer, pointer :: LabRotTable(:,:,:) => null(), LabRotIndx(:) => null ()
integer, allocatable :: vs(:), lrvs(:)
write(*,'("Calculating derivative structures for index n=",i2," to ",i2)') nMin, nMax
write(*,'("Volume",7x,"CPU",5x,"#HNFs",3x,"#SNFs",&
          &4x,"#reduced",4x,"% dups",6x,"volTot",6x,"RunTot")')

! File for timings info
open(99,file="timings_enum.out")

! Set up the output file and write the lattice information
open(14,file="struct_enum.out")
write(14,'(a10)') title
do i = 1,3
   write(14,'(3(g14.8,1x),3x,"# a",i1," parent lattice vector")') parLV(:,i),i
enddo
write(14,'(i2,"-nary case")') k
write(14,'(2i4," # Starting and ending cell sizes for search")') nMin, nMax
if (full) then; write(14,'("Full list of labelings (including incomplete labelings) is used")')
else; write(14,'("Partial list of labelings (complete labelings only) is used")'); endif
write(14,'(8x,"#tot",5x,"#size",1x,"nAt",2x,"pg",4x,"SNF",13x,"HNF",17x,"Left transform",17x,"labeling")')
write(99,'("index",5x,"HNFs",6x,"DupLats",3x,"SNFs",8x,"RemoveDups")')

! Check for 2D or 3D request
if (pLatTyp=='s') then; LatDim = 2
   if (.not. equal(parLV(:,1),(/1._dp,0._dp,0._dp/),eps)) &
        stop 'For "surf" setting, first vector must be 1,0,0'
   if (.not. equal((/parLV(2,1),parLV(3,1)/),(/0._dp,0._dp/),eps)) &
        stop 'For "surf" setting, first component of second and third &
               & must be zero'
else if(pLatTyp=='b') then; LatDim = 3
else; stop 'Specify "surf" or "bulk" in call to "generate_derivative_structures"';endif

! This part generates all the derivative structures and writes the results to the file
runTot = 0
ctot = 0
do ivol = nMin, nMax !max(k,nMin),nMax
   iVolTot = 0
   csize = 0
   call cpu_time(tstart)
   if (LatDim==3) then !<<< 2D ? or 3D? >>
      call get_all_HNFs(ivol,HNF)  ! 3D
   else
      call get_all_2D_HNFs(ivol,HNF) ! 2D
   endif
   call cpu_time(tHNFs)
   call remove_duplicate_lattices(HNF,LatDim,parLV,reducedHNF,fixOp,uqlatts,eps)
   call cpu_time(tDupLat)
   call get_SNF(reducedHNF,L,SNF,B,uqSNF,SNF_labels)
   call cpu_time(tSNF)

   ! We have a list containing each unique derivative superlattice and
   ! its corresponding Smith Normal Form. Now make the labelings.
   Nq = size(uqSNF,3)  ! Number of unique SNFs
   iHNF = 0
   do iuq = 1, Nq ! << Loop over all of the unique SNFs >>
      call cpu_time(tiuq)
      diag = (/uqSNF(1,1,iuq),uqSNF(2,2,iuq),uqSNF(3,3,iuq)/)
      call make_member_list(diag,G)  ! Need the image group G for removing lab-rot dups
      call cpu_time(tML)
      ! Removes trans-dups,non-prims,label-exchange dups
      call generate_labelings(k,diag,labelings,table,trgroup,full) 
      call cpu_time(tGenLab)

      ! Lets create a vector subscript for reducedHNF, L, and fixOp that matches the HNFS
      ! that all have the same SNF
      nRedHNF = size(reducedHNF,3)
      nHNF = count(SNF_labels==iuq)
      if(allocated(vs)) deallocate(vs)
      if(allocated(lrvs)) deallocate(lrvs)
      allocate(vs(nHNF),lrvs(nHNF))
      vs = pack((/(i,i=1,nRedHNF)/),SNF_labels==iuq)
!      write(*,'(20i3)') vs(1:nHNF)
      allocate(LabRotIndx(nRedHNF))
      LabRotIndx = 0
!      SNFmask = reshape(spread(SNF_labels==iuq,1,9),(/3,3,nRedHNF/))
!      HNF = reshape(pack(reducedHNF,SNFmask),(/3,3,nHNF/))
      allocate(LabRotTable(ivol,48,20)) ! Need to find a better way for this...
!      print *,"nHNFs",nHNF
!      do i = 1,nHNF
!         write(*,'(3(3i1,1x),i5)') reducedHNF(:,:,vs(i)), size(fixOp(vs(i))%rot,3)
!      enddo
      call make_label_rotation_table(reducedHNF(:,:,vs),L(:,:,vs),parLV,fixOp(vs),&
                                     G,diag,eps,LabRotTable,LabRotIndx)
!      print *,iuq,"Exited lr_table maker"
      ! Third dimension of LabRotTable is the list of 
!      do i = 1, nHNF
!         do ihnf = 1, ivol
!            write(*,'(8i1)') pack(LabRotTable(:,:,:),LabRotTable(:,:,:)/=0)
!         enddo
!      enddo
!      write(*,'("index: ",20i1)') LabRotIndx(1:nHNF)

!enddo
!enddo
! ******************** Loop over HNF with same perm group
      do ilr = 1, maxval(LabRotIndx) ! loop over the number of label rotation subgroups
         if (associated(tlab)) deallocate(tlab)
         allocate(tlab(size(labelings,1),size(labelings,2)),STAT=status)
         if(status/=0) stop "Allocation of tlab failed in module deriv..."
         tlab = labelings
!         print *, "ilr, max", ilr, maxval(LabRotIndx)
         nrg = count(LabRotIndx==ilr)
!         print *, "nrg",nrg
         call remove_label_rotation_dups(LabRotTable(:,:,ilr),tlab,table,trgroup,k,diag,eps)
         call cpu_time(tRotDup)
         lrvs(1:nHNF) = 0
!         write(*,'("Before pack:", 20i3)') lrvs

         lrvs = pack((/(i,i=1,nHNF)/), LabRotIndx==ilr)
!         write(*,'("Pack:       ", 20i3)') pack((/(i,i=1,nHNF)/),LabRotIndx==ilr)
!         write(*,'("Mask:       ",20l3)') LabRotIndx==ilr
!         write(*,'("LR indx:    ",20i3)') LabRotIndx
!         write(*,'("LR vec sub: ", 20i3)') lrvs
!         write(*,'("vs ",20i3)') vs
         do irg = 1, nrg
            iHNF = iHNF + 1
            ivolTot = ivolTot + size(tlab,1)
            !print *,"iHNF",iHNF
            i = vs(lrvs(irg))
            ld = (/reducedHNF(1,1,i),reducedHNF(2,1,i),reducedHNF(2,2,i),&
                reducedHNF(3,1,i),reducedHNF(3,2,i),reducedHNF(3,3,i)/)
            do ilab = 1,size(tlab,1) ! write out the labelings to a file
               ctot = ctot + 1
               csize = csize + 1
               !is fixOp reference correct here? Yeah, I think so.
               write(14,'(i11,1x,i9,1x,i3,2x,i3,2x,3(i2,1x),2x,6(i2,1x),2x,9(i4),2x,40i1)') &
                    ctot, csize,ivol,size(fixOp(i)%rot,3),diag,ld,&
                    transpose(L(:,:,i)),tlab(ilab,:)
            enddo
         enddo
      enddo ! End of loop over label rotation groups
   enddo
   call cpu_time(tend)
   write(99,'(i2,5x,4(f8.3,2x))') ivol, tHNFs - tstart , tDupLat - tHNFs, tSNF - tDupLat, tend - tGenLab

   
   runTot = runTot + iVolTot
   write(*,'(i4,1x,f14.4,1x,i8,3x,i3,3x,i7,7x,f7.4,i12,i12)')ivol,tend-tstart,size(HNF,3),&
        size(uqSNF,3),size(reducedHNF,3),1-size(reducedHNF,3)/real(size(HNF,3)),ivolTot, runTot
enddo
close(14)
close(99)
END SUBROUTINE generate_derivative_structures

END MODULE derivative_structure_generator

