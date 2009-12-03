MODULE labeling_related
use num_types
use enumeration_types
use vector_matrix_utilities
use numerical_utilities
use rational_mathematics, only: gcd
use combinatorics
implicit none
private
public  get_permutations, count_full_colorings, &
        make_member_list, make_label_rotation_table, generate_unique_labelings, &
        write_labelings
CONTAINS

!***************************************************************************************************
! This routine takes in a list of labels of the parent cell and and the number of different labels
! allowed on each site. It returns a "multiplier", and *expanded* versions of parLabel and parDigit
SUBROUTINE setup_mixed_radix_multiplier(n,k,parLabel,parDigit,label,digit,multiplier)
integer, intent(in) :: n,k ! The index (volume factor) of the supercell, k-nary case
integer, intent(in) :: parLabel(:,:), parDigit(:)
integer, pointer:: label(:,:), digit(:) ! n*longer than parLab/parDigit (INTENT(OUT))
integer, pointer:: multiplier(:)
integer i,j, nD, istat, ic

nD = size(parDigit) ! Size of the d-set n*nD is the total length of a labeling
if (associated(label)) deallocate(label,digit,multiplier)
allocate(label(k,n*nD),digit(n*nD),multiplier(n*nD),STAT=istat)
if (istat/=0) stop "Allocation failed in setup_mixed_radix_multiplier"

digit = (/((parDigit(j),i=1,n),j=1,nD)/) ! Repeat the digit ordinals across all places in the labeling
forall(j=1:k);label(j,:) = (/((parLabel(j,i),ic=1,n),i=1,nD)/); endforall ! Ditto for labels
multiplier = 0; multiplier(n*nD)=1
do i = n*nD-1,1,-1
   multiplier(i) = digit(i+1)*multiplier(i+1)
enddo
END SUBROUTINE setup_mixed_radix_multiplier

!***************************************************************************************************
! This routine takes a list of HNFs, and index, and an indexed list to output the labels for the
! HNFs matching the index. The input for the labelings is not a list of labelings but just the base
! 10 index for those that are unique. The base-10 index is converted to the base-k index and written
! in the output file. Re-expanding the base-10 index to base-k when it was already done once in the
! generate_unique_labelings routine is not efficient CPU-wise but save lots of memory since the
! labelings are never stored in memory except as a base-10 number.
SUBROUTINE write_labelings(k,n,nD,parLabel,parDigit,HNFi,HNFlist,SNFlist,L,fixOp, &
                           Tcnt,Scnt,Hcnt,permIndx,lm,equivalencies,number_ElementN,number_Range)
integer, intent(in) :: k ! number of colors/labels
integer, intent(in) :: n, nD ! index (size of supercell), size of d-set
integer, intent(in) :: parLabel(:,:) ! The *labels* (index 1) for each d-vector (index 2) in the parent
integer, intent(in) :: parDigit(:) ! The *number* of labels allowed on each site of the parent cell 

integer, intent(in) :: HNFi ! Index in the permIndx corresponding to the current block of HNFs
integer, intent(in), dimension(:,:,:) :: HNFlist, SNFlist, L ! List of the HNFs, SNFs, L's. Need this just for the output
type(opList), intent(in) :: fixOp(:)
integer, intent(inout) :: Tcnt, Scnt, Hcnt ! counters for total number of labelings and number of this size and HNF running total 
integer, intent(in) :: permIndx(:)
character, intent(in) :: lm(:)
integer, intent(in) :: number_ElementN
integer, intent(in) :: number_Range(2)

integer nHNF ! Number of HNFs in the list that match the current block index, HNFi
integer, allocatable :: vsH(:), vsL(:) ! Vector subscript for matching the HNF index to the HNF list
integer nl ! Number of unique labelings
integer quot ! quotient for converting base-10 index to base-k labeling
integer iHNF, jHNF, il, i, ilab  ! counters
integer :: labeling(n*nD) ! base-k, n-digit number representing the labeling
integer(li) :: labIndx ! base-10 form of the labeling
integer status ! Allocation exit flag
integer ivsL
integer, pointer :: label(:,:), digit(:), multiplier(:) ! Need to convert base-10 back to labeling

character(3) :: dummy
character(80) :: struct_enum_out_formatstring 


! Labeling Postprocessing data
integer, intent(in)  :: equivalencies(:) ! { full-dset member }: equivalent sites
integer              :: nAllD            ! size of full dset
integer, allocatable :: allD(:)          ! help array: (/ 1, 2, 3, ..., nAllD /)
integer, allocatable :: pplabeling(:)    ! postprocessed labeling
logical :: postprocessLabeling           ! do we need to postprocess labeling
integer, allocatable :: AllDPos2LabelPos(:)  ! { full-dset member }: respective position in labeling

! Postprocessing labelings: setup
nAllD = size(equivalencies)
allocate(allD(nAllD)); allD = (/ (i,i=1,nAllD) /)
allocate(allDPos2LabelPos(nAllD));
! 1) Check whether we have to postprocess the labeling before writing it out
postprocessLabeling = .not. (all(  abs( equivalencies-allD ) ==0))
allocate(pplabeling(n*nAllD))  ! this is ok whether we do postprocessing (nAllD>nD) or not (nAllD==nD)
! 2) if yes:
!    - prepare the postprocessed labeling (pplabeling)
!    - generate a map:  full dset -> enum dset
!      that tells me for a dset member where I can find the correct labeling for it
if (postprocessLabeling) then
  allDPos2LabelPos = (/(count(pack(allD,allD<=equivalencies(i))==equivalencies),i=1,nAllD)/)
  ! this construction is best explained by an example:
  ! suppose we have a dset 1,2,3,4,5
  ! and the equivalent list is 1,4,4,4,5 (so, dset member 1,4,5 will be enumerated, having positions 1,2,3 in labelings)
  ! the map should accomplish: 1 -> 1, 2 -> 2, 3 -> 2, 4 -> 2, 5 -> 3
  ! first, the pack operation take only those in the dset <= equivalent point,
  ! then, the number of truly enumerated points is counted (truly enumerated points are points for which dset=equivalencies
endif



if (any(lm=='')) stop "Labeling index has unmarked entries"
nl = count(lm=='U'); nHNF = count(permIndx==HNFi)
allocate(vsH(nHNF),vsL(nl),STAT=status); if (status/=0) stop "Allocation failed in write_labelings: vsH, vsL"

! Packing...
vsH = pack((/(i,i=1,size(HNFlist,3))/), HNFi==permIndx); 
ivsL=0; do i=1,size(lm); if (lm(i)=='U') then; ivsL=ivsL+1; vsL(ivsL)=i; endif; enddo

! set up the multiplier, labels, digits, etc
call setup_mixed_radix_multiplier(n,k,parLabel,parDigit,label,digit,multiplier)

write(dummy,'(I3)') n*nAllD
struct_enum_out_formatstring = '(i11,1x,i7,1x,i11,1x,i3,2x,i3,2x,3(i2,1x),2x,6(i2,1x),2x,9(i4),2x,'//trim(dummy)//'i1)'
do il = 1, nl ! Loop over the unique labelings
   labIndx = vsL(il)-1 ! Get the base-10 index of the next unique labeling from the vector subscript array
   ! Now convert the base-10 number (labIndx) to the correct labeling
   do ilab=1,n*nD
     quot = labIndx/multiplier(ilab) ! How many times does k(i) divide the number
     labeling(ilab) = label(quot+1,ilab) ! The number of times, indicates the label number
     labIndx = labIndx - quot*multiplier(ilab) ! Take the remainder for the next step
   enddo
   if (postprocessLabeling) then
     call postprocess_labeling(n,nAllD,labeling,pplabeling,allDPos2LabelPos) 
   else
     pplabeling = labeling ! nothing changes
   endif

   !!!do ilab = 1, n*nD ! Loop over each place (digit) in the labeling
   !!!   quot = labIndx/k ! divide the index by k to get the quotient
   !!!   labeling(n*nD-ilab+1) = labIndx - quot*k ! store the remainder (base-k digit of current place)
   !!!   labIndx = quot ! reduce current index to just the quotient and keep going
   !!!enddo
   do iHNF = 1, nHNF ! Write this labeling for each corresponding HNF
      jHNF = vsH(iHNF) ! Index of matching HNFs
      if (check_labeling_numbers(pplabeling,number_ElementN,number_Range)) then
        Tcnt = Tcnt + 1; Scnt = Scnt + 1
        write(14,struct_enum_out_formatstring) &
             Tcnt, Hcnt+iHNF,Scnt,n,size(fixOp(jHNF)%rot,3),SNFlist(1,1,jHNF),SNFlist(2,2,jHNF),SNFlist(3,3,jHNF),&
             HNFlist(1,1,jHNF),HNFlist(2,1,jHNF),HNFlist(2,2,jHNF),HNFlist(3,1,jHNF),HNFlist(3,2,jHNF),&
             HNFlist(3,3,jHNF),transpose(L(:,:,jHNF)),pplabeling   
      endif
   enddo ! loop over HNFs
enddo ! loop over labelings
Hcnt = Hcnt + nHNF

contains

logical function check_labeling_numbers(labeling,el,range)
integer, intent(in) :: labeling(:)
integer, intent(in) :: el
integer, intent(in) :: range(2)
integer             :: c

c = count(labeling==el-1)
if (c>=range(1) .and. c<= range(2)) then
  check_labeling_numbers = .true.
else
  check_labeling_numbers = .false.
endif
end function check_labeling_numbers

subroutine postprocess_labeling(nUC,nAllD,oldlabeling,newlabeling,oldlabel_pos)
integer, intent(in) :: nUC, nAllD
integer, intent(in) :: oldlabeling(:)
integer, intent(out)   :: newlabeling(:)
integer, intent(in) :: oldlabel_pos(:)


integer :: newlab_pos, oldlab_pos
integer :: iD,iUC

do iD=1,nAllD
  do iUC=1,nUC
    newlab_pos = (iD-1)*nUC + iUC
    oldlab_pos = oldlabel_pos(iD)+iUC-1
    newlabeling(newlab_pos) = oldlabeling(oldlab_pos)
  enddo
enddo

end subroutine postprocess_labeling

ENDSUBROUTINE write_labelings

!***************************************************************************************************
! This routine takes in a list of HNFs and a list of rotations (orthogonal transformations) and
! computes the permutations of the labels effected by the rotation. Then the HNFs are grouped into
! categories according to which permutations are applicable.
SUBROUTINE make_label_rotation_table(HNF,L,A,R,G,d, eps, lrTab,lrIndx)
integer, intent(in) :: HNF(:,:,:), L(:,:,:) ! HNFs, set of rotations, left SNF transformations
type(opList), intent(in) :: R(:)
real(dp), intent(in) :: A(3,3) ! Lattice vectors of the parent lattice
integer, intent(in) :: G(:,:) ! The translation group for this SNF
integer, intent(out) :: lrTab(:,:,:) ! Table of the label rotation permutations (perm #, label, lrlist#)
integer, intent(out) :: lrIndx(:) ! Index for permutations list associated with each HNF
integer, intent(in) :: d(3) ! Diagonal elements of the SNF

integer iH, iR, iq, i, j, k, ilq, iM, iHindx, il ! Loop counters
integer nH, nR, n, nlq, nq, status, b, nM(1)
logical unique, err ! flag for identifying new permutation lists, error flag
integer, allocatable :: tlr(:,:) ! temporary list of label rotations for current HNF
real(dp), dimension(3,3) :: T, A1, A1inv, Ainv ! Matrices for making the transformation
integer :: Gp(size(G,1),size(G,2))
integer, allocatable :: trivPerm(:), tM(:,:,:)
real(dp) eps
integer, dimension(3,3) :: M

nH = size(HNF,3) ! Number of HNFs
! Find the maximum number of symmetries for the list of HNFs
nR = 48 ! debug
n = determinant(HNF(:,:,1)) ! Index of the current superlattices
nlq = 0 ! Number of permutation lists that are unique
allocate(tlr(nH,nR),trivPerm(n),STAT=status)
if (status/=0) stop "Trouble allocating tlr, tlrq, trivPerm in make_label_rotation_table"
trivPerm = (/(i,i=1,n)/); tlr = 0;
allocate(tM(3,3,96),STAT=status) ! Why the 96? What's the appropriate number here? >48 but what?
if (status/=0) stop "Trouble allocating tM in make_label_rotation_table"
lrTab = 0; tM = 0; lrIndx = 0

nq = 0; ilq = 0 ! Number of unique transformation matrices (M's), number of unique lists of M's
call matrix_inverse(A,Ainv,err)  ! Need A^-1 to form the transformation
do iH = 1, nH   ! Make a list of permutations for this HNF
   ! First find the permutation on the group by each rotation
   A1 = matmul(L(:,:,iH),Ainv)
   call matrix_inverse(A1,A1inv,err)
   iHindx = 0 ! Keep track of the number of transformations that apply to each HNF
   do iR = 1, size(R(iH)%rot,3) ! Loop over all rotations
      T = matmul(A1, matmul(R(iH)%rot(:,:,iR),A1inv))  ! This is the transformation
      if (.not. equal(T,nint(T),eps)) &
         stop 'ERROR: make_label_rotation_table: Transformation is not integer'
      M = modulo(nint(T),spread(d,2,3))
      ! See if this transformation matrix M is unique
      unique = .true.
      do iq = 1, nq
         if (all(tM(:,:,iq) == M)) then; unique = .false.; exit; endif
      enddo
      if (unique) then
         nq = nq + 1 ! Update number of unique matrices found
         tM(:,:,nq) = M ! Store new matrix in temp list of M's
         iq = nq ! Which matrix in the list is this one? Store in iq
      endif
      if (.not. any(tlr(iH,:)==iq)) then ! This M is not yet associated with this HNF
         iHindx = iHindx + 1             ! so store it in the next place in the list
         tlr(iH,iHindx) = iq
      endif
   enddo ! Loop over rotations
   ! Check if the list of M's for this HNF is unique, or if it is already in the list of lists
   ! First let's sort the list so that the comparisons are robust. (Insertion sort, fine for short lists)
   sort: do j = 2, count(tlr(iH,:)/=0) ! Loop over each non-zero element in the list
      b = tlr(iH,j) ! Temp storage
      k = j-1 ! index pointer to preceding element
      do while(k>0)
         if (tlr(iH,k)<=b) exit
         tlr(iH,k+1)=tlr(iH,k)
         k = k-1
      enddo
      tlr(iH,k+1)=b
   enddo sort ! <<< End sorting
   ! Fail safe on sorting
   do j = 2, count(tlr(iH,:)/=0); if (tlr(iH,j)<tlr(iH,j-1)) stop "Sorting failed";enddo

   ! Now check to see if we have a unique list of M's
   unique = .true.
   do iM = 1, iH
      if (all(tlr(iH,:)==tlr(iM,:)) .and. lrIndx(iM)/= 0) then ! the list has already been found for another HNF
         lrIndx(iH) = lrIndx(iM); unique = .false.; exit; endif
   enddo
   if (unique) then; ilq = ilq + 1; lrIndx(iH) = ilq; endif ! ilq is the number of unique *lists*
enddo

! <<< Make the permutation lists >>>
! We now know which group of permutations is applicable to each HNF. So we just make the lists of
! permutations (right now we just have a list of matrices) and pass that back out
do il = 1, ilq ! Loop over all lists
   ! What is the first list in tlr that corresponds to il?
   nM = minloc(lrIndx,(lrIndx==il))
   do iM = 1, count(tlr(nM(1),:)/=0)
      Gp = matmul(tM(:,:,tlr(nM(1),iM)),G)
      do i=1,3; Gp(i,:) = modulo(Gp(i,:),d(i));enddo  ! Can you do this without the loop, using vector notation?
      do i = 1, n ! Loop over each element of Gp and find its corresponding element in G
         do j = 1, n
            if (all(Gp(:,j)==G(:,i))) then ! the two images are the same for the i-th member
               lrTab(i,iM,il) = j
            endif
         enddo
      enddo
      if (any(lrTab(1:n,iM,il)==0)) then
         stop "Transform didn't work. Gp is not a permutation of G"
      endif
   enddo
enddo

ENDSUBROUTINE make_label_rotation_table

!***************************************************************************************************
! This routine takes two lists of integer sequences and compares them to see if they are the
! same. The input lists are assumed to contain unique entries and perhaps be padded with zeros.
FUNCTION lists_match(list1, list2)
integer, intent(in) :: list1(:,:), list2(:,:) 
logical lists_match, rowmatch
integer nL ! length of each list
integer i,j
lists_match = .false.; rowmatch = .false.
nL = size(list1,2) 
do i = 1, nL
   if (all(list1==list2)) then; rowmatch = .true.; exit; endif
   rowmatch = .false.
   do j = 1, nL
      if (all(list1(:,i)==list2(:,j))) then; rowmatch = .true.; exit; endif
   enddo
   if (.not. rowmatch) exit
enddo
if (rowmatch) lists_match = .true.
ENDFUNCTION lists_match

!***************************************************************************************************
! This function compares a proposed labeling to the list of allowed labels on each site. If any
! label is present "illegally" on one site, then the function is false. This routine is needed
! because in different labels are allowed on different sites, the symmetry operations of the
! underlying parent lattice (where all sites are considered equivalent) can permute a legal labeling
! to an illegal one by permuting one label (on an allowed site) to another site where it is not
! allowed. (GLWH see moleskine 10/9/2009)
FUNCTION labeling_is_legal(labeling,siteLabels,digitN)
logical             :: labeling_is_legal
integer, intent(in) :: labeling(:), siteLabels(:,:), digitN(:)
integer iL, nL

nL = size(labeling)
labeling_is_legal = .true.
!labeling_is_legal = .false.
do iL = 1, nL
   !match = .false.
   if (all(labeling(iL)/=siteLabels(:digitN(iL),iL))) then ! no matching label for this digit
      labeling_is_legal = .false.
      exit ! if there isn't a match in just one digit, illegal labeling
   endif
enddo
!if (match) labeling_is_legal = .true.

ENDFUNCTION labeling_is_legal

!****************************************************************************************************
! This routine takes in the permutations effected by both translation and by rotations that fix the
! superlattice and generates all labelings that are unique. It also removes super-periodic labelings
! (non-primitive superstructures). If the "full" variable is false, it also removes
! "label-permutation" duplicates---labelings that are not unique when the labels themselves (not
! their positions) are permuted (e.g., 00111 <--> 11000).  The basic idea of the routine is to run
! like an "odometer", generating all numbers (base k) from 0 to k^n - 1, and  then use rotation and
! translation permutations to eliminate labelings that represent equivalent superstructures.

SUBROUTINE generate_unique_labelings(k,n,nD,perm,full,lab,parLabel,parDigit)
integer, intent(in) :: k ! Number of colors/labels
integer, intent(in) :: n ! Index of the superlattice
integer, intent(in) :: nD ! Number of sites in the basis of the parent lattice (size of d-set)
integer, intent(in) :: perm(:,:) ! list of translation and rotation permutations
character, pointer :: lab(:) ! Array to store markers for every raw labeling
! I=>incomplete labeling, U=>unique, D=>rot/trans duplicate, N=>non-primitive, E=>label exchange
! Need to pass lab out to write out the labelings
logical, intent(in) :: full ! specify whether the full labelings list should be used or not
integer, intent(in) :: parLabel(:,:) ! The *labels* (index 1) for each d-vector (index 2) in the parent
integer, intent(in) :: parDigit(:) ! The *number* of labels allowed on each site of the parent cell 

integer j ! Index variable (place index) for the k-ary counter
integer ic, i, q ! loop counters, index variables
integer nexp ! number of raw labelings that the k-ary counter should generate
integer nl ! number of labels in each labeling (i.e., determinant size*d-set size)
integer(li) idx ! the base 10 equivalent of the current base k labeling
integer a(n*nD), b(n*nD) ! the "odometer"; label-permuted odometer
integer il ! loop counter over each digit in a labeling
integer digCnt(n*nD) ! Ordinal counter for each place in the labeling (a mixed-radix number)
integer digit(n*nD) ! Each entry is the number of labels in each place
integer label(k,n*nD) ! Same as parLabel but repeated n times
integer multiplier(n*nD) ! place values for each digit. k^(i-1) for the i-th digit for "normal"
!  base-k numbers. More complicated for mixed radix case.
integer c(0:k-1) ! running sum (count) of the number of each label type 
integer id, iq ! Counter for labels that are duplicates, for those unique
integer, pointer :: labPerms(:,:) ! List of permutations of the k labels
integer :: np, ip, nPerm, status ! Loops over label exchang permutations, number of labeling permutatations, allocate error flag

nl = n*nD
nexp = k**nl  ! Number of digits in k-ary counter; upper limit of k-ary counter

!!< Set up the number of expected labelings
nexp = product(parDigit)**n  ! should be the same as k**nl when all labels are on all sites
if (associated(lab)) deallocate(lab)
allocate(lab(nexp),STAT=status)
if(status/=0) stop "Failed to allocate memory for 'lab' in generate_unique_labelings"

! Initialize the counter and ordinal arrays for the mixed-radix counter
digit = (/((parDigit(j),i=1,n),j=1,nD)/) ! Repeat the digit ordinals across all places in the labeling
digCnt = 1 ! Initialize each place to the first label ("lowest" digit)
forall(j=1:k);label(j,:) = (/((parLabel(j,i),ic=1,n),i=1,nD)/); endforall
forall(j=0:k-1); c(j) = count(label(1,:)==j); endforall
!!write(*,'("count",20i2)') c
!!do j = 1,k
!!   write(*,'("initialize label:",1x,20i2)') label(j,:)
!!enddo
!!print *,"nl",nl,"nD",nD,"n",n
!!print *,shape(label)
!!!write(*,'("nl",20i2)') nl
!!!write(*,'("nD",20i2)') nD
!!!write(*,'("mod",20i2)') (mod(i,nD)+1,i=0,nl-1)
!!!
!!write(*,'("digit",20i2)') digit
!!!write(*,'("digCnt",20i2)') digCnt
!!!stop
!!>

a = 0; multiplier = k**(/(i,i=nl-1,0,-1)/) ! The counter; multiplier to convert to base 10

!!< Set up a new multiplier
multiplier = 0; multiplier(nl)=1
!!a = label(1,:); write(*,'("labeling",20i2)') a
do i = nl-1,1,-1
   multiplier(i) = digit(i+1)*multiplier(i+1)
enddo
!!write(*,'("Multiplier: ",10(1x,i6))') multiplier
!!>

lab = ''; iq = 0  ! Index for labelings; number of unique labelings
if (k>12) stop "Too many labels in 'generate_unique_labelings'"
nPerm = size(perm,1)

np = factorial(k) ! Number of permutations of labels (not labelings)

call get_permutations((/(i,i=0,k-1)/),labPerms) 
!call count_full_colorings(k,d,cnt,full) ! Use the Polya polynomial to count the labelings
!call make_translation_group(d,trgrp) ! Find equivalent translations (permutations of labelings)

ic = 0; c = 0; c(0) = nl ! Loop counter for fail safe; initialize digit counter
do; ic = ic + 1
   if (ic > nexp) exit ! Fail safe
   idx = sum(a*multiplier)+1;  ! Index of labeling in base 10
   !!write(*,'(/"index: ",i6)') idx
   !if (idx/=ic) stop "index bug!" !! This check isn't useful for mixed-radix labeling
   if (any(c==0)) then ! Check to see if there are missing digits
      id = id + 1; ! Keep track of the number of incomplete labelings
      if (lab(idx)=='' .and. .not. full) then ! If it isn't marked and we want a partial list, mark it as "incomplete"
         lab(idx) = 'I';    ! Could mark its brothers too...
      endif
   endif
   if (any(c<0)) stop "Bug in the digit counter in generate_unique_labelings: negative digit count!"
   ! If this label hasn't been marked yet, mark it as unique, 'U'
   ! and mark its duplicates as 'D'
   if (lab(idx)=='') then
      lab(idx) = 'U'
      ! Is the first permutation in the list guaranteed to be the identity? We need to skip the identity
      do q = 2,nPerm ! Mark duplicates and eliminate superperiodic (non-primitive) colorings
         if (.not. labeling_is_legal(a(perm(q,:)),label,digit)) then
            !write(*, &
!'("1Skipping illegal labeling",/,"Permutation #: ",i4,/,"original labeling: ",20i1)') q,a
!write(*,'("permuted labeling: ",20i1)') a(perm(q,:))
            !read(*,*)
            cycle
         endif
         idx = sum(a(perm(q,:))*multiplier)+1
         if (idx==ic .and. q <= n) lab(idx)='N' ! This will happen if the coloring is superperiodic
         ! (i.e., non-primitive superstructure). The q=<n condition makes sure we are considering a
         ! "translation" permutation and not a rotation permutation (they're ordered in the
         ! list...and there are n. The first is the identity so skip that one.)
         if (idx > nexp) then
            print *, "Index of permuted labeling is outside the range"
            write(*,'("original labeling ",20i1)') a
            write(*,'("permuted labeling ",20i1)') a(perm(q,:))
            write(*,'("perm ",i2,":",4i2)') (i,perm(i,:),i=1,nPerm)
            stop
         endif
         if (lab(idx)=='') lab(idx) = 'D'  ! Mark as a duplicate
      enddo
      if (.not. full) then ! loop over the label-exchange duplicates and mark them off.
         do q = 1,nPerm ! Loop over all possible permutations. (should this loop start at 1? I think so...)
            do ip = 2,np ! Loop over all permutations of the labels (stored in 'labPerms'). Start at
               ! 2 since we want to skip the identity.
               do il = 1, nl ! For each digit in the labeling, permute the label (label exchange)
                  b(il) = labPerms(a(il)+1,ip)
               enddo
               if (.not. labeling_is_legal(b,label,digit)) then
                  !write(*,'("Permutation #: ",i4,/,"original labeling: ",20i1)') q,a
                  !write(*,'("permuted labeling: ",20i1,"#")') b
                  !read(*,*)
                  cycle
               endif
               !!forall(ia=1:k) ! Convert the k-nary labeling to one with the labels permuted
               !!   where(a==ia-1); b(:) = labPerms(ia,ip);endwhere ! b is the permuted labeling
               !!end forall
               !!if (any(b/=ct)) stop "label exchange duplicate checker is not working..."
               !write(*,'("labeling permuted: ",20i1)') b(perm(q,:))
               idx = sum(b(perm(q,:))*multiplier)+1
               if  (lab(idx)=='') lab(idx) = 'E' ! Only marks of a label-exchange duplicate if it's
               ! not otherwise marked
            enddo
         enddo
      endif ! end block to remove label exchange duplicates
   endif

! "c" counts the number of labels of each kind across the entire labeling. Need this for "partial"
! lists that have label-exchange duplicates removed.
! "digCnt" is the ordinal counter of each digit (i.e., place) in the mixed-radix number (labeling)
   ! Advance the base-k, n*nD-digit counter and keep track of the # of each digit
   j = nl ! Reset the digit index (start all the way to the right again)
   do ! Check to see if we need to roll over any digits, start at the right
      !!if (a(j) /= k - 1) exit ! This digit not ready to roll over, exit the loop and advance digit
      if (digCnt(j) /= digit(j)) exit ! This digit not ready to roll over, exit the loop and advance digit
      !!a(j) = 0  ! Rolling over so set to zero
      a(j) = label(1,j) ! Reset the j-th place to the first digit
      digCnt(j) = 1 ! Reset the ordinal digit count for the j-th place to one
      !!c(k-1) = c(k-1) - 1  ! Update the digit count; just lost a digit of type k-1 (largest possible)
      ! label(digit(j),j) returns the highest number (last symbol) in the j-th place
      c(label(digit(j),j)) = c(label(digit(j),j)) -1 ! Reduce the count of digits of that type
      !!!c(a(j)) = c(a(j)) - 1 ! Reduce the count of digits of that type
      !!c(0) = c(0) + 1 ! So we pick up another zero in the current place ([k-1]->0)
      c(label(1,j)) = c(label(1,j)) + 1  ! So we pick up another "zero" in the current place 
      j = j - 1;       ! Look at the next (going leftways) digit
      if (j < 1) exit  ! If we updated the leftmost digit then we're done
   enddo
   if (j < 1) exit ! We're done counting (hit all possible numbers), exit
   !!a(j) = a(j) + 1 ! Update the next digit (add one to it)
   digCnt(j) = digCnt(j) + 1
   a(j) = label(digCnt(j),j)
   !write(*,'("labeling",20i2)') a
   
   c(a(j)) = c(a(j)) + 1     ! Add 1 to the number of digits of the j+1-th kind
   !!! This doesn't work because the labels aren't necessarily in numerical order
   !!!c(a(j)-1) = c(a(j)-1) - 1 ! subtract 1 from the number of digits of the j-th kind
   c(label(digCnt(j)-1,j)) = c(label(digCnt(j)-1,j)) - 1 ! subtract 1 from the number of digits of the j-th kind
   !write(*,'("iteration:",i10)') ic
   !write(*,'("count",20i2)') c
   !write(*,'("digCnt",20i2)') digCnt
   !write(*,'("j",1x,i1)') j
   if (sum(c) /= nl .and. .not. full) stop 'counting bug'
   !if (ic > 5) stop "early exit for debugging"
enddo
!!!write(*,'(i3,1x,a1)') (i,lab(i),i=1,nexp)
!!!print *,size(lab)
if (ic /= nexp) stop 'Bug: Found the wrong number of labels!'
if (any(lab=="")) stop "Not every labeling was marked in generate_unique_labelings"
END SUBROUTINE generate_unique_labelings

 
!**************************************************************************************************
! Takes the length of three cyclic group and constructs the member list so that
! the permutations can be determined. Each member has three components, corresponding to the
! entries for each of the three cyclic groups. 
SUBROUTINE make_member_list(n,p)
integer, intent(in)  :: n(3)! Diagonal elements of the SNF
integer, pointer :: p(:,:)  ! List of members of the translation group

integer im, status  ! loop over members, allocate status flag
if (associated(p)) deallocate(p)
allocate(p(3,product(n)),STAT=status)
if(status/=0) stop "Allocation of 'p' failed in make_member_list"
p = 0
do im = 2,product(n)  ! Loop over the members of the translation group
   p(:,im) = p(:,im-1) ! Start with the same digits as in the previous increment
   p(3,im) = mod(p(3,im-1)+1,n(3))  ! Increment the first cyclic group
   if (p(3,im)==0) then             ! If it rolled over then 
      p(2,im) = mod(p(2,im-1)+1,n(2))! increment the next cyclic group
      if (p(2,im)==0) then          ! If this one rolled over too
         p(1,im) = mod(p(1,im-1)+1,n(1)) ! Then increment the third one
      endif
   endif
enddo
ENDSUBROUTINE make_member_list

!***************************************************************************************************
! This routine finds all the permutations of the group members that leave the decoration unchanged.
! Essentially we are finding a list of mappings: add to the group one of the members of the group
! to get another member of the group. This list of mappings is the list of labeling (coloring) 
! permutations that leave the superstructure unchanged.
SUBROUTINE make_translation_group(d,trans)
integer, intent(in) :: d(3) ! members of the group, diagonal elements of SNF
integer, pointer :: trans(:,:) ! Translations that leave the superstructure unchanged

integer n, im,i, j, status
integer tg(3,d(1)*d(2)*d(3)) ! Temporary storage for the translations
integer, pointer :: m(:,:) => null() ! List of group members

! Get a list of group members; essentially a mixed-radix, 3 digit counter
call make_member_list(d,m)
n = product(d) ! Number of elements in each permutation group
if (associated(trans)) deallocate(trans)
allocate(trans(n,n),STAT=status)
if(status/=0) stop "Allocation of 'trans' failed in make_translation_group"

do im = 1, n
   ! add the im-th element of the group to every element in the group, mod d
   forall(i=1:size(tg,2)); tg(:,i) = mod(m(:,im) + m(:,i),d);end forall
   ! Find the index of the group member in the translated group
   do i = 1,n ! This approach is an N^2 loop. Can this be improved? Does it matter? I don't think so N is small
      do j = 1,n
         if (all(tg(:,i)==m(:,j))) then ! the two members are equal
            trans(im,i) = j             ! Save the index and exit the loop
            exit; endif; enddo; enddo
   if (sum(trans(im,:)) /= (n*n+n)/2) stop "*** Bug in group generator! ***"
enddo
END SUBROUTINE make_translation_group

!*******************************************************************************
! This routine calculates the result of the counting polynomial to find the total number
! of unique colorings of "n" things using "k" colors. The outer loop is for inclusion/exclusion
! so that colorings that don't use each color at least once are removed from the counting. The use
! of every color at least once is what is meant by a "full" coloring
SUBROUTINE count_full_colorings(k,d,count,full)
integer, intent(in) :: k, d(3)  ! Number of colors, diagonal entries of the Smith Normal Form matrix
integer, intent(out) :: count   ! Number of unique, full colorings
logical, intent(in) :: full ! Use a full list? (including incomplete labelings)
integer x1,x2,x3 ! Counters over d1,d2,d3, in the triple sum
integer p        ! Counter over number of terms in the exclusion/inclusion counting
integer m, tc    ! Color index (counter), intermediate (temporary) counter
integer n        ! Number of elements to be colored
m = k; count = 0
n = product(d)
if (full) then   ! The number of labelings for a full list is trivial, just k^n
   count = k**n
else  ! For a partial list (incomplete labelings not included) things are slightly more complicated
   do p = 0,m-1  ! Loop over each term in the exclusion/inclusion expression
      tc = 0
      do x1 = 1,d(1)  ! Loop over the triple sum
      do x2 = 1,d(2)
      do x3 = 1,d(3)
         tc = tc + (m-p)**gcd(d(1)*d(2)*d(3), x1*d(2)*d(3), d(1)*x2*d(3), d(1)*d(2)*x3)
      enddo; enddo; enddo
      count = count + (-1)**p*nchoosek(m,p)*tc
   enddo
   count = count/n
endif

ENDSUBROUTINE count_full_colorings
END MODULE labeling_related
