! Generate all derivative structures of a parent structure
! Gus L. W. Hart BYU July 2007
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
use sorting
implicit none
private
public get_all_HNFs, remove_duplicate_lattices, get_SNF, get_all_2D_HNFs,&
     &  gen_multilattice_derivatives, &
     get_dvector_permutations, get_rotation_perms_lists, do_rotperms_form_groups
CONTAINS

!***************************************************************************************************
! This function checks that every "product" of two permutations in a list is still in the list. That
! is that every list forms a group. I'm not sure this is a sufficient condition but it is a
! necessary one. (Better to check not for closure only (as here) but to check that for each member
  ! of the group A there is the inverse of another element, B^-1, whose composition AB^-1 leads to a
  ! third element C that is a member of the supposed group.)
FUNCTION do_rotperms_form_groups(rpl)
logical do_rotperms_form_groups
type(RotPermList) :: rpl(:) ! The reduced (unique) set of rotation permutations lists for a given index

integer nL, iL, iP, jP, nP, kP, ng, itest
logical exists
integer, allocatable :: testperm(:)

nL = size(rpl) ! Number of lists
ng = size(rpl(1)%perm,2) ! Number of elements in each permutations
allocate(testperm(ng))

do_rotperms_form_groups = .true.
lists: do iL = 1, nL ! Check each list in the set
   nP = size(rpl(iL)%perm,1) ! Number of permutations in this list
   write(11,'(3x)',advance="no")
   do iP = 1, nP; write(11,'(i3)',advance="no") iP; enddo; write(11,*) ! Make the column headings
perms:   do iP = 1, nP
      write(11,'(i3)',advance="no") iP   ! Write the row heading
      !do jP = 1, nP; write(11,'(3x)',advance='no');enddo ! skip to the column
      do jP = 1, nP ! Now loop over all products for the iP'th permutation
         exists = .false. 
         ! Is the product of perm_i x perm_j in the set?
         ! Permute the elements of the iP'th permutation according to the jP'th permutation
         testperm = rpl(iL)%perm(iP,rpl(iL)%perm(jP,:))
         do kP = 1, nP
            if (all(testperm==rpl(iL)%perm(kp,:))) then
               exists = .true.
               write(11,'(i3)',advance="no") kP
               exit
            endif
         enddo
         if (.not. exists) then
            
            write(11,'(/,"The set of permutations doesn''t form a group")')
            write(11,'("failed at ",2i3)')  jP, iP
            write(11,'("testperm ",20i2,/)') testperm
            write(11,'("permlist ",8i2)') transpose(rpl(iL)%perm)
            !write(*,'(20i2,/)') testperm
            print *
            !do itest = 1,nP
            !   write(*,'(20i2)') rpl(iL)%perm(itest,:)
            !enddo
            do_rotperms_form_groups = .false.
            !exit lists
            exit perms
         endif
      enddo
      write(11,*)
   enddo perms
   if (exists)write(11,*) iL,"The permutations have closure"
   if(.not. exists) write(11,*) iL,"failed"
enddo lists ! Loop over lists
!if(exists) then; do_rotperms_form_groups = .true.
!else; do_rotperms_form_groups = .false.
!endif
ENDFUNCTION do_rotperms_form_groups
!***************************************************************************************************
! This routine takes a list of permutation lists and identifies those that are idential. The
! routine assumes that lists associated with different SNFs should be kept separate. The output,
! RPLindx is an index groups the lists that match. The lists themselves are assumed to be in
! "alphabetical" order so that they can be quickly compared.
SUBROUTINE organize_rotperm_lists(SNFindx, RPList,rdRPList,RPLindx)
integer, intent(in):: SNFindx(:)
type(RotPermList), intent(in) :: RPList(:)
type(RotPermList), pointer :: rdRPList(:) ! output
integer, pointer :: RPLindx(:) ! output

integer iL, jL, nL, status, iSNF, nSNF, cnt, Scnt
type(RotPermList), allocatable :: tList(:)
logical unique

nL = size(SNFindx); nSNF = maxval(SNFindx)
allocate(tList(nL),RPLindx(nL),STAT=status)
if(status/=0) stop "Allocation failed in organize_rotperm_lists: tList, RPLindx"

Scnt = 0
do iSNF = 1, nSNF ! Do a separate treatment for each SNF form
   cnt = 0
   do iL = 1, count(SNFindx==iSNF)
      unique = .true.
      do jL = 1, cnt ! Loop over the number of unique lists for this SNF
         ! if the lists aren't the same length, then they definitely aren't identical
         if (size(RPList(iL)%perm,1)/=size(tList(Scnt + jL)%perm,1)) cycle
         if (all(RPList(iL)%perm==tList(Scnt + jL)%perm)) then
            unique = .false.
            exit
         endif
      enddo
      if (unique) then ! store this list in the master list
         cnt = cnt + 1
         nL = size(RPList(iL)%perm,1)
         allocate(tList(cnt + Scnt)%perm(nL,size(RPList(iL)%perm,2)))
         tList(cnt + Scnt)%nL = nL
         tList(cnt + Scnt)%perm = RPlist(iL)%perm
      endif
      ! Store the label for this unique list in the output index
      RPLindx(iL) = jL + Scnt ! Not sure this is going to do the right thing...
   enddo
   Scnt = Scnt + cnt
enddo
! Now copy the reduced list to the output variable
allocate(rdRPList(Scnt))
do iL=1,Scnt
   nL = size(tList(iL)%perm,1)
   allocate(rdRPList(iL)%perm(nl,size(tList(iL)%perm,2)));
   rdRPList(iL)%nL = nL
   rdRPList(iL)%perm = tList(iL)%perm
enddo
ENDSUBROUTINE organize_rotperm_lists

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
!   do iD = 1,3
!      write(*,'(i2,": ",3i2)') iOp,nint(rot(iD,:,iOp))
!   enddo
!   write(*,'("S:",3i3,/)') nint(shift(:,iOp))
   tRD = rd
   do iD = 1, nD
      call bring_into_cell(rd(:,iD),inv_pLV,pLV,eps)
   enddo
   dRPList%v(:,:,iOp) = rd(:,:) - tRD(:,:)
   call map_dvector_permutation(rd,pd,dRPList%perm(iOp,:),eps)
enddo
! I don't think we should reduce this list to a unique one. Some rotations that don't permute the d's could
! still permute the g's. So we have to keep all the d permutations, even if they look redundant here.

ENDSUBROUTINE get_dvector_permutations
!***************************************************************************************************
! For each HNF, we have a list of the operations (rotations + shifts, if present) that leave the
! superlattice fixed. Given this set of fixing operations, make a list of the permutations on the
! d-vectors (interior points of the multilattice) effected by the rotations. Then sort the HNFs into
! blocks that have the same rotation permutations.
SUBROUTINE get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps)
real(dp), intent(in) :: A(3,3) ! Lattice vectors of the primary lattice (parent lattice)
integer, intent(in), dimension(:,:,:) :: HNF, L, SNF ! List of HNF matrices, left transforms, and their SNFs
type(OpList), intent(in) :: Op(:) ! A list of symmetry ops (rots and shifts) for the parent multilattice
type(RotPermList) :: RPlist(:) ! A list of lists of permutations effected by the Ops
type(RotPermList), intent(in) :: dperms
real(dp), intent(in) :: eps ! finite precision tolerance
type(RotPermList) :: rperms
integer, pointer :: g(:,:)
integer, allocatable :: gp(:,:), dgp(:,:) ! G prime; the "rotated" group, (d',g') "rotated"  table
integer, allocatable :: tg(:,:), perm(:), ident(:,:) ! translated group, translation permutation of the group members
integer iH, nH, diag(3), iD, nD, iOp, nOp, n, ig, it, status, im, jm
real(dp), dimension(3,3) :: Ainv, T, Tinv
logical err
real(dp), allocatable :: rgp(:,:)
logical, allocatable :: skip(:)

nH = size(HNF,3); n = determinant(HNF(:,:,1)); nD = size(RPList(1)%v,2) ! Num superlattices, index, Num d's
allocate(gp(3,n), dgp(nD,n), rgp(3,n),skip(n),STAT=status)
if(status/=0) stop "Allocation failed in get_rotation_perm_lists: gp, dgp, rgp, skip" 
allocate(tg(3,n),perm(n),ident(nD,n),STAT=status)
if(status/=0) stop "Allocation failed in get_rotation_perm_lists: tg, perm, ident" 
ident = transpose(reshape((/(ig,ig=1,n*nD)/),(/n,nD/)))
forall(iH = 1:nH); RPlist(iH)%nL=0; end forall

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
   Tinv = matmul(L(:,:,iH),Ainv); call matrix_inverse(Tinv, T, err)
   if (err) stop "Bad inverse for transformation matrix: get_rotation_perm_lists"
   nOp = size(Op(iH)%rot,3);
   write(*,'(3i2)') transpose(HNF(:,:,iH))
   if (associated(rperms%perm)) deallocate(rperms%perm)
   allocate(rperms%perm(nOp*n+n+nOp,nD*n),STAT=status)
   if (status/=0) stop "Allocation failed in get_rotation_perm_lists: rperms%perm"
   do iOp = 1, nOp
      dgp = 0
      do iD = 1, nD ! Loop over each row in the (d,g) table
         ! LA^-1(v_i+(RAL^-1)G)
         rgp = matmul(Tinv,(spread(RPList(iH)%v(:,iD,iOp),2,n)+matmul(matmul(Op(iH)%rot(:,:,iOp),T),g)))
         if (.not. equal(rgp,nint(rgp),eps)) stop "Transform left big fractional parts"
         gp = nint(rgp)
         gp = modulo(gp,spread(diag,2,n)) ! Mod by each entry of the SNF to bring into group
         ! Now that the permuted group has been found, find the mapping of the elements between the
         ! orginal group and the permuted group. This is the permutation.
         
!         call find_permutation_of_group(g,gp,dgp(dperms%perm(RPList(iH)%RotIndx(iOp),iD),:))
!         dgp(dperms%perm(RPList(iH)%RotIndx(iOp),iD),:) = dperms%perm(RPList(iH)%RotIndx(iOp),iD)+(iD-1)*n
         skip = .false. ! This is just for efficiency
         do im = 1, n
            do jm = 1, n
               if (skip(jm)) cycle
               if (all(gp(:,jm)==g(:,im))) then
                  dgp(dperms%perm(RPList(iH)%RotIndx(iOp),iD),im) = jm+(iD-1)*n
                  skip(jm) = .true.
                  exit
               endif
            enddo ! jm
         enddo ! im
      enddo ! loop over d-vectors
      if (any(dgp==0)) stop "(d,g)-->(d',g') mapping failed in get_rotation_perm_lists"

      ! Now we have the (d',g') table for this rotation. Now record the permutation
      rperms%perm(iOp,:) = reshape(transpose(dgp),(/nD*n/)) ! "fast index" is iD, then g
      rperms%nL = nOp
   enddo ! loop over rotations 

   ! Now that we have the permutations that are effected by N+t type of rotations (for each N in the
   ! point group), we need to compose them with the permutations effected by lattice translations,
   ! (parent) lattice translations contained inside the supercell (N+t+r, in Rod's
   ! nomenclature). Only when we have these two lists, and compose them, do we have a guarantee that
   ! the resulting list is actually a group.

   ! To get the permutations effected by the r's, we don't need the r's. We can merely take the
   ! member list, the group, and add each element of the group to the group itself and see what
   ! permutation happens. (This part is somewhat redundant with make_translation_group in
   ! labeling_related module.)
   do ig = 1,n ! The number of r's inside the superlattice (the translation perms) is the same as the index n
      tg = g(:,:)+spread(g(:,ig),2,n)
      tg = modulo(tg,spread(diag,2,n)) ! mod by the SNF entries
      call find_permutation_of_group(g,tg,perm)
!print *,shape(spread(perm,2,nD))
!print *,shape(spread(perm,1,nD))
      rperms%perm(nOp+ig,:) = reshape(transpose(ident(:,perm)),(/n*nD/))
   ! Now we need to compose the N+t permutations with r permutations (rotations composed with
   ! translations) to get a list that should be group
 
   enddo
   do iOp = 1,nOp ! Loop over rotation perms (N+t type)
      do it = 1,n ! Loop overe translation perms (r type)
         ! Form the permutation effected by composing the iOp-th one with the it-th one
         rperms%perm(nOp+n+(iOp-1)*n+it,:) = rperms%perm(iOp,(rperms%perm(nOp+it,:)))
         write(*,'("const",30i2)') rperms%perm(nOp+n+iOp+it-1,:)
      enddo
   enddo

   print *,"total list"
   do ig = 1,nOp+n+nOp*n
      write(*,'("tot",i5,1x,30i3)') ig,rperms%perm(ig,:)
   enddo
!   stop
   ! Sort the permutations and add them to the master list
   call sort_permutations_list(rperms%perm)
   ! The permutations list is now in "alphabetical" order and contains no duplicates, so allocate
   ! the next list in RPList and store this one. The duplicate lists, if any, will be identified in
   ! a separate step
      print *,"sorted list"
   do ig = 1,size(rperms%perm,1)
      write(*,'(30i3)') rperms%perm(ig,:)
   enddo
   print *,"Number of ops",nOp
   print *," n",n
   print *,'iH',iH
   !if (iH==3) stop
!stop
   RPlist(iH)%nL = size(rperms%perm,1)
   allocate(RPlist(iH)%perm(RPlist(iH)%nL,n*nD))
   RPlist(iH)%perm = rperms%perm
enddo ! loop over iH (superlattices)
ENDSUBROUTINE get_rotation_perms_lists

!***************************************************************************************************
SUBROUTINE find_permutation_of_group(g,gp,perm)
integer, intent(in), dimension(:,:) :: g, gp ! unpermuted and permuted groups
integer, intent(out) :: perm(:) ! permutation of gp

integer n ! number of elements in the group (index of the superlattice)
logical skip(size(g,2))
integer im, jm
n = size(g,2); perm = 0
skip = .false. ! This is just for efficiency
do im = 1, n
   do jm = 1, n
      if (skip(jm)) cycle ! This is just for efficiency
      if (all(gp(:,jm)==g(:,im))) then
         perm(im) = jm
         skip(jm) = .true.
         exit ! don't keep looking if you already found the match
      endif
   enddo ! jm
enddo ! im
!write (*,'(30i3)') perm
if (any(perm==0)) stop "mapping failed in find_permutation_of_group"
ENDSUBROUTINE find_permutation_of_group

!***************************************************************************************************
SUBROUTINE map_dvector_permutation(rd,d,RP,eps)
real(dp), dimension(:,:) :: rd, d
integer :: RP(:)
real(dp), intent(in) :: eps

integer iD, jD, nD
logical found(size(RP))

RP = 0; found = .false.
nD = size(rd,2) ! # of d-vectors
do iD = 1, nD
   do jD = 1, nD
      if(found(jD)) cycle
      if(equal(rd(:,iD),d(:,jD),eps)) then
         RP(iD)=jD
         found(jD) = .true.

         exit
      endif
   enddo
enddo
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
SUBROUTINE get_SNF(HNF,A,SNF,B,RPList,F,SNF_label,fixing_op)
integer, pointer :: HNF(:,:,:), SNF_label(:)
integer, pointer, dimension(:,:,:) :: A, SNF, B, F
type(opList) :: fixing_op(:) ! List of operations that fix each HNF
type(RotPermList) :: RPList(:) ! List of rotation permutations for each HNF

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
RPList = RPList(indx)

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
SUBROUTINE remove_duplicate_lattices(hnf,LatDim,parent_lattice,d,dperms,uq_hnf,fixing_op,RPList,latts,eps)
integer, pointer :: hnf(:,:,:) ! HNF matrices (input)
integer :: LatDim ! Is the parent lattice 2D or 3D?
real(dp), intent(in) :: parent_lattice(3,3) ! parent lattice (input)
integer, pointer :: uq_hnf(:,:,:) ! list of symmetrically distinct HNFs (output)
type(opList), pointer :: fixing_op(:) ! List of operations that fix each HNF
real(dp),intent(in):: eps ! finite precision (input)
real(dp), pointer :: d(:,:) 
type(RotPermList), intent(in) :: dperms
type(RotPermList), pointer :: RPList(:)

real(dp), pointer:: sgrots(:,:,:), sgshift(:,:)
real(dp), dimension(3,3) :: test_latticei, test_latticej, thisRot, rotLat, origLat
integer i, Nhnf, iuq, irot, j, nRot, Nq, ic, status, nD
integer, allocatable :: temp_hnf(:,:,:), tIndex(:)
real(dp), pointer :: latts(:,:,:)
logical duplicate
type(opList), allocatable :: tmpOp(:)
real(dp), allocatable :: tSGrots(:,:,:), tv(:,:,:)
integer, pointer :: aTyp(:)


nD = size(d,2)
allocate(aTyp(nD),STAT=status)
if(status/=0)stop "Allocation failed in remove_duplicate_lattices: aTyp"
aTyp = 1

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
allocate(uq_hnf(3,3,iuq),latts(3,3,iuq),RPList(iuq),STAT=status)
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
allocate(tv(3,nD,nRot),tIndex(nRot),STAT=status); if(status/=0) stop "tv didn't allocate"
do i=1,Nq; allocate(tmpOp(i)%rot(3,3,nRot),tmpOp(i)%shift(3,nRot),STAT=status)
if(status/=0) stop "Allocation failed in remove_duplicat_lattices: tmpOp%rot or shift";enddo
do iuq = 1, Nq  ! Loop over each unique HNF
   ic = 0
   tv = 0; tIndex = 0
   do iRot = 1,nRot  ! Loop over each rotation
      thisRot = sgrots(:,:,iRot) ! Store the rotation
      origLat = matmul(parent_lattice,uq_hnf(:,:,iuq))  ! Compute the superlattice
      rotLat = matmul(thisRot,origLat)          ! Compute the rotated superlattice
      if (is_equiv_lattice(rotLat,origLat,eps)) then ! this operation fixes the lattice and should be recorded
         ic = ic + 1
         tmpOp(iuq)%rot(:,:,ic) = thisRot
         tmpOp(iuq)%shift(:,ic) = sgshift(:,iRot)
         tv(:,:,ic) = dperms%v(:,:,iRot)
         tIndex(ic) = iRot
      endif
   enddo ! Now we know which rotations fix the lattice and how many there are so store them
   do i=1,ic; allocate(fixing_op(iuq)%rot(3,3,ic),fixing_op(iuq)%shift(3,ic),STAT=status)
      if(status/=0) stop "Allocation of fixing_op(iuq) failed, module deriv..."; enddo ! Allocate the storage for them
   fixing_op(iuq)%rot =   tmpOp(iuq)%rot(:,:,1:ic) ! Stuff the rotations into the permanent array
   fixing_op(iuq)%shift = tmpOp(iuq)%shift(:,1:ic) ! Stuff the shifts into the permanent array
   allocate(RPList(iuq)%v(3,nD,ic),RPList(iuq)%RotIndx(ic),STAT=status)
   if (status/=0) stop "Allocation of RPList failed in remove_duplicate_lattices" 
   RPList(iuq)%v = tv(:,:,1:ic)
   RPList(iuq)%RotIndx = tIndex(1:ic)
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
integer, pointer, dimension(:,:,:) :: HNF => null(),SNF => null(), L => null(), R => null()
integer, pointer :: labelings(:,:) =>null(), SNF_labels(:) =>null(), tlab(:,:)=>null(), uqSNF(:,:,:) => null()
integer, pointer, dimension(:,:,:) :: rdHNF =>null()
real(dp) tstart, tend
type(opList), pointer :: fixOp(:)  ! Symmetry operations that leave a multilattice unchanged
type(RotPermList), pointer :: RPList(:) ! Master list of the rotation permutation lists
type(RotPermList), pointer :: rdRPList(:) ! Master list of the *unique* rotation permutation lists
type(RotPermList) :: ParRPList ! Just the list for the parent lattice
integer, pointer ::  RPLindx(:) => null() ! Index showing which list of rotation permutations corresponds to which HNF
real(dp), pointer :: uqlatts(:,:,:) => null()

!debug variables
real(dp), dimension(3,3) :: sLV, sLVi ! superlattice and inverse
real(dp), dimension(:,:), pointer :: fulld=>null(),rd=>null() !dset, rotated
integer, pointer :: rp(:,:)=>null() ! rotation perm
integer itest, jd
logical err
type(rotpermlist), pointer :: testlist(:)
!end debug variables

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
call get_dvector_permutations(parLV,d,ParRPList,eps)
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
   call remove_duplicate_lattices(HNF,LatDim,parLV,d,ParRPList,rdHNF,fixOp,RPList,uqlatts,eps)
   ! Superlattices with the same SNF will have the same list of translation permutations of the
   ! labelings. So they can all be done at once if we find the SNF. rdHNF is the reduced list.
   call get_SNF(rdHNF,L,SNF,R,RPList,uqSNF,SNF_labels,fixOp)
   ! Each HNF will have a certain number of rotations that leave the superlattice fixed, called
   ! fixOp. These operations will effect a permutation on the (d,g) table. Since many of the HNFs
   ! will have an identical list of rotation permutations, it'll be efficient to reduce the
   ! labelings for all such HNFs just once. So we need to generate the list for each HNF and then
   ! sort the HNFs into blocks with matching permutation lists.
   call get_rotation_perms_lists(parLV,rdHNF,L,SNF,fixOp,RPList,ParRPList,eps)
   call organize_rotperm_lists(SNF_labels,RPList,rdRPList,RPLindx)

!   write(*,'("RPLindx ",30i3)') RPLindx
!   print *,rplindx
!   write(*,'(20i2)') transpose(rdHNF(:,:,3))
!   itest = 3
!   sLV = matmul(parLV,HNF(:,:,itest))
!   call matrix_inverse(sLV,sLVi,err)
!   if (err) stop "matrix inverse failed"
!   !do iD = 1, 3
!   allocate(fulld(3,nD*ivol),rd(3,nD*ivol),rp(16,nd*ivol))
!  
!  call expand_dvectors(ivol,parLV,sLV,sLVi,HNF(:,:,itest),SNF(:,:,itest),L(:,:,itest),d,fulld,eps)
!   ! write out the dvectors
! do jd = 1, nD*iVol
!    write(*,'(3i3)') nint(fulld(:,jd))
! enddo;print *
!   do iD = 1, size(fixOp(itest)%rot,3) ! loop over rotations
!      ! rotate fulld
!      rd = matmul(fixOp(itest)%rot(:,:,iD),fulld)+spread(fixOp(itest)%shift(:,iD),2,nD)
!  
!      write(*,'("Rot no.: ",i3)') iD
!      do jd = 1, size(fulld,2)
!         call bring_into_cell(rd(:,jd),sLVi,sLV,eps)
!      enddo
!      if (iD ==13) then
!         write(*,'(3i3)') nint(transpose(fixOp(itest)%rot(:,:,iD))),nint(fixOp(itest)%shift(:,iD)) 
!         
!    do jd = 1, nD*iVol
!       write(*,'(3i3)') nint(rd(:,jd))
!    enddo; print *
!endif

      ! find mapping
!      call map_dvector_permutation(rd,fulld,rp(iD,:),eps)
!  
!      write(*,'(20i2)') rp(iD,:)
!   enddo
!   call sort_permutations_list(rp); print *
!   write(*,'("master ",8i2)') transpose(rp(:,:))
!   allocate(testlist(1))
!   allocate(testlist(1)%perm(8,8))
!   testlist(1)%perm = rp
!   print *, do_rotperms_form_groups(testlist)
!   stop "end of test"
   !write(*,'(20i2)') nint(transpose(fixOp(3)%rot(:,:,4)))
   !write(*,'(20i2)') nint((fixOp(3)%shift(:,4)))
   !
   !write(*,'(20i2)') nint(transpose(fixOp(3)%rot(:,:,6)))
   !write(*,'(20i2)') nint((fixOp(3)%shift(:,6)))

   if (.not. do_rotperms_form_groups(rdRPList)) print *, "Rotperm list doesn't form group"
enddo ! loop over cell sizes (ivol)
call cpu_time(tend)

ENDSUBROUTINE gen_multilattice_derivatives
END MODULE derivative_structure_generator

