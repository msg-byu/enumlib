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
     & generate_derivative_structures, gen_multilattice_derivatives, &
     get_dvector_permutations
CONTAINS

!***************************************************************************************************
SUBROUTINE get_dvector_permutations(pLV,pd,dRPList,eps)
real(dp) :: pLV(3,3) ! Lattice vectors of the primary lattice (parent lattice)
real(dp), pointer :: pd(:,:) ! d-vectors defining the multilattice (primary lattice only)
type(RotPermList), intent(out) :: dRPList ! Output. A list permutations effected by the Ops
real(dp), intent(in) :: eps ! finite precision tolerance

integer nD, iD, nOp, iOp, status
integer, pointer :: aTyp(:), tList(:,:)
real(dp) :: rd(size(pd,1),size(pd,2)), tRD(size(pd,1),size(pd,2))
real(dp) :: inv_pLV(3,3) ! Inverse of the pLV matrix
real(dp), pointer:: rot(:,:,:), shift(:,:), tv(:,:,:)
logical err

nD = size(pd,2)
allocate(aTyp(nD),STAT=status)
if(status/=0)stop "Allocation failed in get_dvector_permutations: aTyp"
aTyp = 1

call get_spaceGroup(pLV,aTyp,pd,rot,shift,.false.,eps)
nOp = size(rot,3)
allocate(tList(nOp,nD),tv(3,nD,nOp),STAT=status)
if(status/=0)stop "Allocation failed in get_dvector_permutations: tList"
allocate(dRPList%perm(nOp,nD),dRPList%v(3,nD,nOp),STAT=status)
if(status/=0)stop "Allocation failed in get_dvector_permutations: dRPList"


call matrix_inverse(pLV,inv_pLV,err)
if (err) stop "Bad parent lattice vectors in input to get_dvector_permutations"

do iOp = 1, nOp ! Try each operation in turn and see how the d-vectors are permuted for each
   rd = matmul(rot(:,:,iOp),pd)+spread(shift(:,iOp),2,nD) ! Rotate each d and add the shift
   tRD = rd
   do iD = 1, nD
      call bring_into_cell(rd(:,iD),inv_pLV,pLV,eps)
   enddo
   dRPList%v(:,:,iOp) = rd(:,:) - tRD(:,:)
   call map_dvector_permutation(rd,pd,dRPList%perm(iOp,:),eps)
enddo

! I don't think making the unique list is correct. Some rotations that don't permute the d's could
!still permute the g's. So we have to keep all the d permutations, even if they look redundant here.

!! Now that we have the permutations, we need to reduce the list to only those that are unique
!! Since this is a small set (no bigger than 48) and we only do this once, a straightforward approach
!! should be fine (that is, efficiency isn't an issue here).
!nqP = 0
!do iOp = 1, nOp
!   unique = .true.
!   do iqP = 1, nqP
!      if (all(dRPlist%perm(iOp,:)==tList(iqP,:))) then
!         unique = .false.
!         exit
!      endif
!   enddo ! End loop over unique permutations
!   if (unique) then ! found a new permutation
!      nqP = nqP + 1
!      tList(nqP,:) = dRPlist%perm(iOp,:)
!      tv(:,:,nqP) = dRPList%v(:,:,iOp)
!   endif
!enddo ! loop over operations
!deallocate(dRPList%perm,dRPList%v)
!allocate(dRPList%perm(nqP,nD),dRPList%v(3,nD,nqP))      
!dRPList%perm = tList(1:nqP,:)
!dRPList%v = tv(:,:,1:nqP)

!do iqP = 1, nqP
!   write(*,'(20i1)') dRPList%perm(iqP,:)
!enddo
!do iqP = 1, nqP
!   write(*,'(20(3i3,1x))') nint(dRPList%v(:,:,iqP))
!   !write(*,'(20(3(f7.3,1x)))') dRPList%v(:,:,iqP)
!enddo

ENDSUBROUTINE get_dvector_permutations
!***************************************************************************************************
! For each HNF, we have a list of the operations (rotations + shifts, if present) that leave the
! superlattice fixed. Given this set of fixing operations, make a list of the permutations on the
! d-vectors (interior points of the multilattice) effected by the rotations. Then sort the HNFs into
! blocks that have the same rotation permutations.
SUBROUTINE get_rotation_perm_lists(A,HNF,L,SNF,dperms,Op,RPlist,RPLx,eps)
real(dp), intent(in) :: A(3,3) ! Lattice vectors of the primary lattice (parent lattice)
integer, intent(in), dimension(:,:,:) :: HNF, L, SNF ! List of HNF matrices, left transforms, and their SNFs
type(OpList), intent(in) :: Op(:) ! A list of symmetry ops (rots and shifts) for the parent multilattice
type(RotPermList), pointer :: RPlist(:) ! Output. A list of lists of permutations effected by the Ops
integer, pointer :: RPLx(:) ! An index indicating which HNF is subject to which list of permutations
real(dp), intent(in) :: eps ! finite precision tolerance
type(RotPermList), intent(in) :: dperms 
type(RotPermList) :: rperms
integer, pointer :: dgPermList(:,:), g(:,:), dg(:,:,:)
integer, allocatable :: gp(:,:), dgp(:,:) ! G prime; the "rotated" group, (d',g') "rotated" table
integer iH, nH, diag(3), iD, nD, iOp, nOp, n, im, jm, status, nL
real(dp), dimension(3,3) :: Ainv, T, Tinv
logical err
logical, allocatable :: skip(:)
real(dp), allocatable :: rgp(:,:)

nH = size(HNF,3); n = determinant(HNF(:,:,1)); nD = size(dperms%v,2) ! Num superlattices, index, Num d's
allocate(gp(3,n), dgp(nD,n), skip(n),rgp(3,n),STAT=status)
if(status/=0) stop "Allocation failed in get_rotation_perm_lists: gp, dgp, skip, rgp" 
allocate(RPLx(nH),RPlist(nH))
forall iH = 1:nH; RPlist(iH)%nL=0; endforall

! Make the group member list for the first SNF
diag = (/SNF(1,1,1),SNF(2,2,1),SNF(3,3,1)/)
call make_member_list(diag,g)
call matrix_inverse(A,Ainv,err)
if(err) stop "Invalid parent lattice vectors in get_rotation_perm_lists"

do iH = 1,nH ! loop over each superlattice
   ! unless the SNF is different than the previous (they should be sorted into blocks) don't bother
   ! making the group again. Just use the same one.
   if(iH > 1) then; if(.not. all(SNF(:,:,ih)==SNF(:,:,ih-1))) then
      diag = (/SNF(1,1,ih),SNF(2,2,ih),SNF(3,3,ih)/)
      call make_member_list(diag,g)
   endif; endif
   !if (all(HNF(:,:,iH)==reshape((/1,1,0,0,2,0,0,0,2/),(/3,3/)))) print *,"########" 
   Tinv = matmul(L(:,:,iH),Ainv); call matrix_inverse(Tinv, T, err)
   if (err) stop "Bad inverse for transformation matrix: get_rotation_perm_lists"
   nOp = size(Op(iH)%rot,3);!print*,"nop",nOp
   if (associated(rperms%perm)) deallocate(rperms%perm)
   allocate(rperms%perm(nOp,nD*n))
   do iOp = 1, nOp
!write(*,'(3(f7.3,1x))') transpose(Op(iH)%rot(:,:,iOp));print*
      dgp = 0
      do iD = 1, nD ! Loop over each row in the (d,g) table
         !write(*,'("iH:",i2,3x,"iOp:",i2,3x,"iD:",i2)') iH,iOp,iD
!write(*,'(3(i2,1x))') transpose( SNF(:,:,iH));print*
!write(*,'(3(f7.3,1x))') transpose( matmul(A,HNF(:,:,iH)));print*
!write(*,'(4(f7.3,1x))') spread(dperms%v(:,iD,iOp),2,n); print *
!write(*,'(4(f7.3,1x))') matmul(Tinv,(spread(dperms%v(:,iD,iOp),2,n)+matmul(matmul(Op(iH)%rot(:,:,iOp),T),g))); print*
!write(*,'(3(f7.3,1x))') transpose(oP(iH)%rot(:,:,iOp));print*
         ! LA^-1(v_i+(RAL^-1)G)
         rgp = matmul(Tinv,(spread(dperms%v(:,iD,iOp),2,n)+matmul(matmul(Op(iH)%rot(:,:,iOp),T),g)))
         if (.not. equal(rgp,nint(rgp),eps)) stop "Transform left big fractional parts"
         gp = nint(rgp)
         !write(*,'(4i2)') transpose(gp); print *
         !write(*,'(4i2)') transpose(spread(diag,2,n));print*
         gp = modulo(gp,spread(diag,2,n)) ! Mod by each entry of the SNF to bring into group
!write(*,'(4i2)') transpose(gp); print*
         skip = .false.
         do im = 1, n
            do jm = 1, n
               if (skip(jm)) cycle
               if (all(gp(:,jm)==g(:,im))) then
 !                 write(*,'("iOp",i3,1x,"loc:",4i3)')iop,dperms%perm(iOp,iD)
                  dgp(dperms%perm(iOp,iD),im) = jm+(iD-1)*n
                  skip(jm) = .true.
                  exit
               endif
            enddo ! jm
         enddo ! im
  !       write(*,'(4(i2,1x))') transpose(dgp); print *
      enddo ! loop over d-vectors
      !write(*,'(4(i2,1x))') transpose(dgp-1)
      if (any(dgp==0)) stop "(d,g)-->(d',g') mapping failed in get_rotation_perm_lists"
      ! Now we have the (d',g') table for this rotation. Now record the permutation
      rperms%perm(iOp,:) = reshape(transpose(dgp),(/nD*n/))
      rperms%nL = nOp
      !write(*,'(16i3)') rperms%perm(iOp,:)-1
   enddo ! loop over rotations
   ! Now see if this list of rotation permutations is unique. Put it in the master list if so.
   call add_perms_to_master_list(rperms,RPlist,PRLx,nL)
enddo ! loop over iH (superlattices)

! Get rid of the empty entries is RPlist

ENDSUBROUTINE get_rotation_perm_lists
!***************************************************************************************************
! This routine takes a list of permutations and compares it to a list of other lists of
! permutations. If it is unique, it is added to the master list.
SUBROUTINE add_perms_to_masterlist(inList, Master, indx, nL)
type(RotPermList), intent(in) :: inList, Master(:)
integer, pointer, intent(out) :: indx(:)
integer, intent(out) :: nL ! Number of permutation lists in the master list

integer iL, nL

call sort_perms_list(inList)
nL = size(master)
unique = .true.
do iL = 1, nL
   if (master%nL==0) exit

enddo
if(unique) then ! put this in the list
   master(iL) = in
endif
ENDSUBROUTINE add_perms_to_masterlist


!!***************************************************************************************************
!! For each HNF, we have a list of the operations (rotations + shifts, if present) that leave the
!! superlattice fixed. Given this set of fixing operations, make a list of the permutations on the
!! d-vectors (interior points of the multilattice) effected by the rotations. Then sort the HNFs into
!! blocks that have the same rotation permutations.
!SUBROUTINE get_rotation_perm_lists(pLV,pd,HNF,L,SNF,Op,RPlist,RPLx,eps)
!real(dp) :: pLV(3,3) ! Lattice vectors of the primary lattice (parent lattice)
!real(dp), pointer :: pd(:,:) ! d-vectors defining the multilattice (primary lattice only)
!integer, intent(in), dimension(:,:,:) :: HNF, L, SNF ! List of HNF matrices, left transforms, and their SNFs
!type(OpList), intent(in) :: Op(:) ! A list of symmetry ops (rots and shifts) for the parent multilattice
!type(RotPermList), pointer :: RPlist(:) ! Output. A list of lists of permutations effected by the Ops
!integer, pointer :: RPLx(:) ! An index indicating which HNF is subject to which list of permutations
!real(dp), intent(in) :: eps ! finite precision tolerance
!
!integer iR, iH, iD, idummy ! various counters
!integer n, nD ! Index of superlattice, number of d-vectors
!integer status ! Flag for allocations that fail
!real(dp), pointer :: d(:,:), rd(:,:), trd(:,:) ! d-vectors of the superlattice, rotated versions of the same
!real(dp) :: sLV(3,3), sLVinv(3,3)
!logical err
!integer, pointer :: tRP(:), tRPlist(:,:) ! temp. rotation perm, list of the same
!real(dp), pointer :: v(:,:,:)  !vector that returns the rotated d-vector back into the first cell
!
!n = determinant(HNF(:,:,1)) ! Index of the superlattices in this block
!nD = n*size(pd,2) ! 
!allocate(d(3,nD),Rd(3,nD),trd(3,nD),STAT=status)
!if(status/=0) stop "Allocation failed in get_rotation_perm_lists: d, rd, trd"
!allocate(tRP(nD),tRPlist(48,nD), v(3,nD,48),STAT=status)
!if(status/=0) stop "Allocation failed in get_rotation_perm_lists: tRP, tRPlist, v"
!
!! For each HNF, loop over all fixing Ops and hit the multilattice with the Op. Then move the
!! d-vectors back to the first unit cell. Then construct an index vector (the permutation) that
!! indicates how each d-vector in the rotated lattice matches with those in the original.
!do iH = 1,size(HNF,3)
!   !print *,"********************************"
!   !write(*,'(/,20(3i3,/))') transpose(HNF(:,:,iH))
!   
!   sLV = matmul(pLV,HNF(:,:,iH))
!   call matrix_inverse(sLV,sLVinv,err)
!   if(err) stop "Matrix inverse didn't work in get_rotation_perm_lists"
!
!   ! Need to define the d-vectors of the superlattice (right now we only have those for the parent)
!   call expand_dvectors(n,pLV,sLV,sLVinv,HNF(:,:,iH),SNF(:,:,iH),L(:,:,iH),pd,d,eps)
!   do iR = 1,size(Op(iH)%rot,3)
!      idummy = size(Op(iH)%rot,3)
!      !write(*,'("nRot:",i3," dummy")') idummy
!      !write(10,'("nRot:",i3)') size(Op(iH)%rot,3)
!      !print *,"dummy"
!
!      ! Rotate all the d-vectors by one of the rotations and add the shift, if present
!      rd = matmul(Op(iH)%rot(:,:,iR),d)+spread(Op(iH)%shift(:,iR),2,nD)
!      trd = rd ! Make a temporary copy before we bounce the rd's back into the cell
!      !write(*,'(/,20(3i3,/))') nint(d(:,1:nD))
!      !write(*,'(/,20(3i3,/))') nint(rd(:,1:nD))
!      !
!      !write(*,'(/,20(3f7.3,/))') Op(iH)%shift(:,iR)
!      !write(*,'(/,20(3f7.3,/))') matmul(sLVinv,rd(:,1:nD))
!      do iD = 1, nD ! Move each of the d-vectors into the first unit cell of the superlattice
!         call bring_into_cell(rd(:,iD),sLVinv,sLV,eps)
!      enddo 
!      v(:,:,iR) = rd(:,:) - trd(:,:)
!      !write(*,'(/,"Rot#: ",i2,/)') iR
!      !write(*,'(/,20(3f7.3,/))') matmul(sLVinv,rd(:,1:nD))
!
!      !write(*,'(/,20(3i3,/))') nint(rd(:,1:nD))
!
!      ! determine how the d-vectors have been permuted by the symmetry operation
!      !write(*,'(/,20(3i3,/))') transpose(nint(Op(iH)%Rot(:,:,iR)))
!      call map_dvector_permutation(rd,d,tRPlist(iR,:),eps)
!   enddo ! loop over rotations
!   !write(*,'("HNF:",3(3i3,/),/)') transpose(HNF(:,:,iH))
!   
!   
!   ! tRPlist now has a list of permutations. We need to remove redundant perms and then
!   ! order the list so it is easy to compare to other lists for uniqueness. Probably ought to remove
!   ! instances of the trivial permutation too, for the sake of efficiency.
!   !call make_permutations_unique(tRPlist,v)
!   !call add_permutations_to_list
!   
!enddo ! end loop over HNFs
!ENDSUBROUTINE get_rotation_perm_lists

!!***************************************************************************************************
!! This routine takes a list of permutations, removes the trivial permutation and all other
!! permutations that are redundant. 
!SUBROUTINE make_permutations_unique(list,v)
!integer, pointer :: list(:,:)
!
!integer nP, nD, iP, jP, iD, trivPerm(size(list,2)), uqP, nq, status
!logical found, unique
!integer, allocatable :: tlist(:,:) ! Temporary copy of "list"
!real(dp), allocatable :: tv(:,:,:) ! Temporary copy of the v's
!
!nP = size(list,1); nD = size(list,2)
!trivPerm = (/(iD,iD=1,nD)/); nq = 0
!allocate(tlist(nP,nD),tv(3,nD,nP),STAT=status)
!if (status/=0) stop "Allocation failed in make_permutations_unique: tlist"
!do iP = 1, nP
!   if (all(list(iP,:) == trivPerm)) cycle
!   unique = .true.
!   do uqP = 1, nq
!      if (all(list(iP,:)==tlist(uqP,:))) then
!         found = .false. ! Found
!         exit
!      endif
!   enddo ! loop over unique permutations
!   if (unique) then
!      nq = nq + 1
!      tlist(nq,:) = list(iP,:)
!      tv(:,:,nq) = v(:,:,iP)
!   endif
!enddo ! loop over permutations
!! store unique permutations
!deallocate(list,v)
!allocate(list(nq,nD),v(3,nD,nq),STAT=status)
!if (status/=0) stop "Allocation failed in make_permutations_unique: list"
!list = tlist(1:nq,:) ! Store the reduced list
!v = tv(:,:,1:nq)
!ENDSUBROUTINE make_permutations_unique

!***************************************************************************************************
SUBROUTINE map_dvector_permutation(rd,d,RP,eps)
real(dp), dimension(:,:) :: rd, d
integer(si) :: RP(:)
real(dp), intent(in) :: eps

integer iD, jD, nD
logical found(size(RP))

RP = 0; found = .false.
nD = size(rd,2) ! # of d-vectors
!write(*,'("Rot",20(3i2,/),/)') nint(rd)
!write(*,'("Org",20(3i2,/),/)') nint(d)
do iD = 1, nD
   do jD = 1, nD
      if(found(jD)) cycle
      !write(*,'(3f7.3,1x,3f7.3,2x,l1)')rd(:,iD),d(:,jD),equal(rd(:,iD),d(:,jD),eps)
      if(equal(rd(:,iD),d(:,jD),eps)) then
         RP(iD)=jD
         found(jD) = .true.
         exit
      endif
   enddo
   !print *
enddo
!write(*,'("perm: ",4i2)') RP-1

if(any(RP==0)) then; print *, "d-vector didn't permute in map_dvector_permutation";stop;endif


ENDSUBROUTINE map_dvector_permutation

!***************************************************************************************************
! Take a superlattice (parent + HNF) and a set of interior points and generate all the multilattice
! points in the supercell
SUBROUTINE expand_dvectors(n,LV,sLV,sLVi,H,S,L,pd,d,eps)
real(dp) :: LV(3,3), eps ! Lattice vectors of the parent cell, finite precision tolerance
integer H(3,3), S(3,3), L(3,3) ! HNF, SNF matrix, and the left transform for the given HNF
integer n ! Index (size) of the supercell
real(dp), pointer :: pd(:,:), d(:,:) ! interior points of the primary lattice, points of the multilattice
real(dp), intent(in), dimension(3,3) :: sLV, sLVi ! superlattice vectors and inverse

integer, pointer :: g(:,:)
integer id, Linv(3,3)
real(dp) Li(3,3)
logical err

!write(*,'(3(3i3,/),/)') transpose(H)
!write(*,'(3(3i3,/),/)') transpose(L)

call make_member_list((/S(1,1),S(2,2),S(3,3)/),g)
call matrix_inverse(real(L,dp),Li,err)
if(.not. equal(Li,nint(Li),eps)) stop "Transform didn't work in expand d-vectors"
Linv = nint(Li)

d(:,1:n) = matmul(LV,matmul(Linv,g)) ! Get the set of superlattice d-vectors that are at the origin of
                                     ! each parent cell in the superlattice 

do id = 0, size(pd,2)-1 
      d(:,id*n+1:id*n+n) = d(1:3,1:n) + spread(pd(1:3,id+1),2,n) 
enddo ! loop over group elements
! d now contains n copies of pd, each copy shifted by the origin of each parent cell in the superlattice
do id = 1, n*size(pd,2) ! Move each of the primary cell d-vectors into the first unit cell
   call bring_into_cell(d(:,id),sLVi,sLV,eps)
enddo ! loop over primary cell d-vectors


!write(*,'(/,20(3i3,/))') nint(d(:,1:n))

ENDSUBROUTINE expand_dvectors

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
aTyp = 1

call get_spaceGroup(parent_lattice,aTyp,d,sgrots,sgshift,.false.,eps)
nRot = size(sgrots,3)
!print *,"nRot", nRot
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
! This routine is should eventually replace "generate_derivative_structures". The difference with
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
integer, pointer, dimension(:,:,:) :: HNF => null(),SNF => null(), L => null(), B => null()
integer, pointer :: labelings(:,:) =>null(), SNF_labels(:) =>null(), tlab(:,:)=>null(), uqSNF(:,:,:) => null()
integer, pointer, dimension(:,:,:) :: rdHNF =>null()
real(dp) tstart, tend
type(opList), pointer :: fixOp(:)  ! Symmetry operations that leave a multilattice unchanged
type(RotPermList), pointer :: RPList(:) ! Master list of the rotation permutation lists
type(RotPermList) :: ParentDvecRotPermList ! Just the list for the parent lattice
integer, pointer ::  RPLindx(:) ! Index showing which list of rotation permutations corresponds to which HNF
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
!write(14,'("Symmetry of the primary lattice is of order ",i2)')
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

! The permutations of the interior points (d-vectors) under symmetry operations of the parent
! multilattice are used later on. Generate them here
call get_dvector_permutations(parLV,d,ParentDvecRotPermList,eps)
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
   call get_SNF(rdHNF,L,SNF,B,uqSNF,SNF_labels,fixOp)
   ! Each HNF will have a certain number of rotations that leave the superlattice fixed, called
   ! fixOp. These operations will effect a permutation on the (d,g) table. Since many of the HNFs
   ! will have an identical list of rotation permutations, it'll be efficient to reduce the
   ! labelings for all such HNFs just once. So we need to generate the list for each HNF and then
   ! sort the HNFs into blocks with matching permutation lists.
   !call get_rotation_perm_lists(parLV,d,rdHNF,L,SNF,fixOp,RPList,RPLindx,eps)
   call get_rotation_perm_lists(parLV,rdHNF,L,SNF,ParentDvecRotPermList,fixOp,RPList,RPLindx,eps)
    !call
enddo ! loop over cell sizes (ivol)
call cpu_time(tend)

ENDSUBROUTINE gen_multilattice_derivatives
END MODULE derivative_structure_generator

