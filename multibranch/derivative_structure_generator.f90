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
public get_all_HNFs, remove_duplicate_lattices, get_SNF, get_all_2D_HNFs,&
     & generate_derivative_structures, gen_multilattice_derivatives
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
! corresponds to each HNF. The transformations, SNFs, and labels are ordered on output
! into blocks of identical SNFs.
SUBROUTINE get_SNF(HNF,A,SNF,B,F,SNF_label,fixing_op)
integer, pointer :: HNF(:,:,:), SNF_label(:)
integer, pointer, dimension(:,:,:) :: A, SNF, B, F
type(opList), pointer :: fixing_op(:) ! List of operations that fix each HNF

integer :: indx(size(HNF,3))
integer ihnf, nHNF, nfound, ifound, status, ic, jc, i
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

! We have all the SNFs now, as well as the transformations. But they are not in any kind of
! order. Rearrange the transformations, labels, and HNFs so that they are in blocks of matching SNFs
ic = 1; jc = 0
do ifound = 1, nfound
   jc = count(SNF_label==ifound) ! How many HNFs have this SNF?
   indx(ic:ic+jc-1) = pack((/(i,i=1,nHNF)/),SNF_label==ifound) ! Get the index of each matching case
   ic = ic + jc ! Compute the offset where the index is stored
enddo

! Use the computed order index, indx, to reorder the A, B, and SNF matrices
HNF = HNF(1:3,1:3,indx)
SNF = SNF(1:3,1:3,indx)
A = A(1:3,1:3,indx)
B = B(1:3,1:3,indx)
SNF_label = SNF_label(indx) ! Reorder the labels for which SNF each HNF is
fixing_op = fixing_op(indx)

! Fail safe trigger and storage of unique SNFs
if (ic/=nHNF+1) stop "SNF sort in get_SNF didn't work"
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
SUBROUTINE remove_duplicate_lattices(hnf,LatDim,parent_lattice,d,uq_hnf,fixing_op,latts,eps)
integer, pointer :: hnf(:,:,:) ! HNF matrices (input)
integer :: LatDim ! Is the parent lattice 2D or 3D?
real(dp), intent(in) :: parent_lattice(3,3) ! parent lattice (input)
integer, pointer :: uq_hnf(:,:,:) ! list of symmetrically distinct HNFs (output)
type(opList), pointer :: fixing_op(:) ! List of operations that fix each HNF
real(dp),intent(in):: eps ! finite precision (input)
real(dp), pointer :: d(:,:)

real(dp), pointer:: sgrots(:,:,:), sgshift(:,:)
real(dp), dimension(3,3) :: test_latticei, test_latticej, thisRot, rotLat, origLat
integer i, Nhnf, iuq, irot, j, nRot, Nq, ic, status, nD
integer, allocatable :: temp_hnf(:,:,:)
real(dp), pointer :: latts(:,:,:)
logical duplicate
type(opList), allocatable :: tmpOp(:)
real(dp), allocatable :: tSGrots(:,:,:)
integer, pointer :: aTyp(:)

nD = size(d,2)
allocate(aTyp(nD),STAT=status)
if(status/=0)stop "Allocation failed in remove_duplicate_lattices: aTyp"

call get_spaceGroup(parent_lattice,aTyp,d,sgrots,sgshift,.false.,eps)
nRot = size(sgrots,3)

Nhnf = size(hnf,3)
allocate(temp_hnf(3,3,Nhnf),STAT=status)
if(status/=0) stop "Failed to allocate memory in remove_duplicate_lattices: temp_hnf"
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
if(status/=0) stop "Allocation of tmpOp or fixing_op failed in remove_duplicate_lattices"
do i=1,Nq; allocate(tmpOp(i)%rot(3,3,nRot),tmpOp(i)%shift(3,nRot),STAT=status)
if(status/=0) stop "Allocation failed in remove_duplicat_lattices: tmpOp%rot or shift";enddo
do iuq = 1, Nq  ! Loop over each unique HNF
   ic = 0
   do iRot = 1,nRot  ! Loop over each rotation
      thisRot = sgrots(:,:,iRot) ! Store the rotation
      origLat = matmul(parent_lattice,uq_hnf(:,:,iuq))  ! Compute the superlattice
      rotLat = matmul(thisRot,origLat)          ! Compute the rotated superlattice
      if (is_equiv_lattice(rotLat,origLat,eps)) then ! this operation fixes the lattice and should be recorded
         ic = ic + 1
         tmpOp(iuq)%rot(:,:,ic) = thisRot
         tmpOp(iuq)%shift(:,ic) = sgshift(:,iRot)
      endif
   enddo ! Now we know which rotations fix the lattice and how many there are so store them
   do i=1,ic; allocate(fixing_op(iuq)%rot(3,3,ic),fixing_op(iuq)%shift(3,ic),STAT=status)
if(status/=0) stop "Allocation of fixing_op(iuq) failed, module deriv..."; enddo ! Allocate the storage for them
   fixing_op(iuq)%rot = tmpOp(iuq)%rot(:,:,1:ic) ! Stuff the rotations into the permanent array
   fixing_op(iuq)%shift = tmpOp(iuq)%shift(:,1:ic) ! Stuff the shifts into the permanent array
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
SUBROUTINE generate_derivative_structures(title, parLV, nD, d, k, nMin, nMax, pLatTyp, eps, full)
integer, intent(in) :: k, nMin, nMax 
character(10), intent(in) :: title
real(dp), intent(in) :: parLV(3,3)
character(1), intent(in) :: pLatTyp
logical, intent(in) :: full 

integer, pointer, dimension(:,:,:) :: uqSNF => null()
integer, pointer :: trgroup(:,:)=>null()
real(dp), pointer :: uqlatts(:,:,:) => null(), d(:,:)=>null()
character, pointer :: table(:)
integer i, ihnf, ilab, iuq, ivol, Nq, iVolTot, runTot, LatDim, status, ilr, nHNF, nRedHNF
integer irg, nrg, nD
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
   call remove_duplicate_lattices(HNF,LatDim,parLV,d,reducedHNF,fixOp,uqlatts,eps)
   call cpu_time(tDupLat)
   call get_SNF(reducedHNF,L,SNF,B,uqSNF,SNF_labels,fixOp)
   call cpu_time(tSNF)
   write(*,'(20i3)') SNF_labels
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

!***************************************************************************************************
! This routine should eventually replace "generate_derivative_structures". The difference with
! this one is that it applies to superstructures derived from multilattices. This is more general
! and should therefore work on "mono"-lattices, as the original routine did. The algorithm
! implemented here has been slightly reorder from the original routine and that discussed in the
! first paper.
SUBROUTINE gen_multilattice_derivatives(title, parLV, nD, d, k, nMin, nMax, pLatTyp, eps, full)
integer, intent(in) :: k, nMin, nMax, nD 
character(10), intent(in) :: title
real(dp), intent(in) :: parLV(3,3), eps
real(dp), pointer :: d(:,:)
character(1), intent(in) :: pLatTyp
logical, intent(in) :: full 

integer iD, i, ivol, LatDim, runTot, ctot, ivoltot, csize
integer, pointer, dimension(:,:,:) :: HNF => null(), rHNF => null(),SNF => null(), L => null(), B => null()
integer, pointer :: labelings(:,:) =>null(), SNF_labels(:) =>null(), tlab(:,:)=>null(), uqSNF(:,:,:) => null()
real(dp) tstart, tend
type(opList), pointer :: fixOp(:)  ! Symmetry operations that leave a multilattice unchanged
type(RPlist), pointer :: rotPermList(:) ! Master list of the rotation permutation lists
integer, pointer ::  RPLindx ! Index showing which list of rotation permutations corresponds to which HNF
real(dp), pointer :: uqlatts(:,:,:) => null()


write(*,'("Calculating derivative structures for index n=",i2," to ",i2)') nMin, nMax
write(*,'("Volume",7x,"CPU",5x,"#HNFs",3x,"#SNFs",&
          &4x,"#reduced",4x,"% dups",6x,"volTot",6x,"RunTot")')

! Set up the output file and write the lattice information
open(14,file="struct_enum.out")
write(14,'(a10)') title
do i = 1,3
   write(14,'(3(g14.8,1x),3x,"# a",i1," parent lattice vector")') parLV(:,i),i
enddo
write(14, '(i5," # Number of points in the multilattice")') nD
do iD = 1,nD
   write(14,'(3(g14.8,1x),3x,"# d",i2.2," d-vector")') d(:,iD),iD
enddo
write(14,'(i2,"-nary case")') k
write(14,'(2i4," # Starting and ending cell sizes for search")') nMin, nMax
if (full) then; write(14,'("Full list of labelings (including incomplete labelings) is used")')
else; write(14,'("Partial list of labelings (complete labelings only) is used")'); endif
write(14,'(8x,"#tot",5x,"#size",1x,"nAt",2x,"pg",4x,"SNF",13x,"HNF",17x,"Left transform",17x,"labeling")')

! Check for 2D or 3D request
if (pLatTyp=='s') then; LatDim = 2
   if (.not. equal(parLV(:,1),(/1._dp,0._dp,0._dp/),eps)) &
        stop 'For "surf" setting, first vector must be 1,0,0'
   if (.not. equal((/parLV(2,1),parLV(3,1)/),(/0._dp,0._dp/),eps)) &
        stop 'For "surf" setting, first component of second and third &
               & must be zero'
else if(pLatTyp=='b') then; LatDim = 3
else; stop 'Specify "surf" or "bulk" in call to "generate_derivative_structures"';endif

! This part generates all the derivative structures. Results are writen to unit 14.
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
   ! Many of the superlattices will be symmetrically equivalent so we use the symmetry of the parent
   ! multilattice to reduce the list to those that are symmetrically distinct.
   call remove_duplicate_lattices(HNF,LatDim,parLV,d,rdHNF,fixOp,uqlatts,eps)
   ! Superlattices with the same SNF will have the same list of translation permutations of the
   ! labelings. So they can all be done at once if we find the SNF. rHNF is the reduced list.
   call get_SNF(rHNF,L,SNF,B,uqSNF,SNF_labels,fixOp)
   ! Each HNF will have a certain number of rotations that leave the superlattice fixed, call
   ! fixOp. These operations will effect a permutation on the (d,g) table. Since many of the HNFs
   ! will have an identical list of rotation permutations, it'll be efficient to reduce the
   ! labelings for all such HNFs just once. So we need to generate the list for each HNF and then
   ! sort the HNFs into blocks with matching permutation lists.
   call get_rotation_perm_lists(rHNF,L,SNF,fixOp,rotPermList,RPLindx)

enddo ! loop over cell sizes (ivol)

ENDSUBROUTINE gen_multilattice_derivatives
END MODULE derivative_structure_generator

