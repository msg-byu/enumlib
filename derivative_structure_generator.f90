! Generate all derivative structures of a parent structure
! Gus L. W. Hart BYU July 2007
! Things still to do:
! * Currently only applies to primitive parent structures*

MODULE derivative_structure_generator
use num_types
use crystal_types
use combinatorics
use symmetry_module
use compare_structures
use rational_mathematics
use numerical_utilities
use vector_matrix_utilities, only: determinant, matrix_inverse

implicit none
private
public get_all_HNFs, remove_duplicate_lattices, get_SNF, get_all_2D_HNFs
CONTAINS

!***************************************************************************************************
! Finds all the possible diagonals of the HNF matrices of a given size
SUBROUTINE get_HNF_diagonals(detS, diagonals)
integer, intent(in) :: detS ! Cell size, i.e., determinant of S matrix
integer, pointer :: diagonals(:,:) ! All possible diagonals

integer i, j, id, quotient
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
allocate(diagonals(3,id))
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
allocate(F(3,3,nfound))
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
   allocate(sgrots(3,3,nRot))
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
allocate(uq_hnf(3,3,iuq),latts(3,3,iuq))
uq_hnf = temp_hnf(:,:,:iuq)
Nq = iuq
forall (i=1:Nq); latts(:,:,i) = matmul(parent_lattice,uq_hnf(:,:,i));end forall
! Now loop through the list of unique superlattices and record which of the parent lattice
! symmetry operations leave the superlattice unchanged. These operations, while they leave
! the superlattice unchanged can permute the labels inside. Thus, this is another source of 
! duplicates. 
allocate(tmpOp(Nq),fixing_op(Nq))
do i=1,Nq; allocate(tmpOp(i)%rot(3,3,nRot));enddo
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
   do i=1,ic; allocate(fixing_op(iuq)%rot(3,3,ic)); enddo ! Allocate the storage for them
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

integer i, id, quotient
integer :: tempDiag(3,detS*3)

id = 0 ! Number of diagonals found
do i = 1,detS ! Loop over possible first factors
   if (.not. mod(detS,i)==0) cycle
   quotient = detS/i
   id = id + 1
   tempDiag(:,id) = (/1, i, quotient /) ! Construct the factor triplet
enddo
allocate(diagonals(3,id))
diagonals = tempDiag(:,1:id)
END SUBROUTINE get_HNF_2D_diagonals

END MODULE derivative_structure_generator

