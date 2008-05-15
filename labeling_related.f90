MODULE labeling_related
use num_types
use enumeration_types
use vector_matrix_utilities
use numerical_utilities
use rational_mathematics, only: gcd
use combinatorics
implicit none
private
public  remove_label_rotation_dups, get_permutations, count_full_colorings, &
       generate_labelings, make_member_list, make_label_rotation_table
CONTAINS
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

integer iH, iR, jR, iq, i, j, k, ilq, iM, iHindx, il ! Loop counters
integer nH, nR, n, nlq, nq, status, tNr, b, nM(1)
logical unique, err ! flag for identifying new permutation lists, error flag
integer, allocatable :: tlr(:,:) ! temporary list of label rotations for current HNF
real(dp), dimension(3,3) :: T, A1, A1inv, Ainv ! Matrices for making the transformation
integer :: Gp(size(G,1),size(G,2))
integer, allocatable :: trivPerm(:), tM(:,:,:)
real(dp) eps
integer, dimension(3,3) :: M

!print *,"shape HNF in: ",shape(HNF)
nH = size(HNF,3) ! Number of HNFs
! Find the maximum number of symmetries for the list of HNFs
!nR = 0; do i = 1, size(R); tNr = size(R(i)%rot,3); if(tNr>nR) nR = tNr; enddo  
nR = 48 ! debug
n = determinant(HNF(:,:,1)) ! Index of the current superlattices
nlq = 0 ! Number of permutation lists that are unique
allocate(tlr(nH,nR),trivPerm(n),STAT=status)
if (status/=0) stop "Trouble allocating tlr, tlrq, trivPerm in make_label_rotation_table"
trivPerm = (/(i,i=1,n)/); tlr = 0;
allocate(tM(3,3,96),STAT=status)
if (status/=0) stop "Trouble allocating tM in make_label_rotation_table"
lrTab = 0; tM = 0; lrIndx = 0

nq = 0; ilq = 0 ! Number of unique transformation matrices (M's), number of unique lists of M's
do iH = 1, nH   ! Make a list of permutations for this HNF
   ! First find the permutation on the group by each rotation
   call matrix_inverse(A,Ainv,err)  ! Need A^-1 to form the transformation
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
   do j = 2, count(tlr(iH,:)/=0) ! Loop over each non-zero element in the list
      b = tlr(iH,j) ! Temp storage
      k = j-1 ! index pointer to preceding element
      do while(k>0)
         if (tlr(iH,k)<=b) exit
         tlr(iH,k+1)=tlr(iH,k)
         k = k-1
      enddo
      tlr(iH,k+1)=b
   enddo ! <<< End sorting
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
      if (any(lrTab(1:n,iM,il)==0)) stop "Transform didn't work. Gp is not a permutation of G"
   enddo
enddo

ENDSUBROUTINE make_label_rotation_table

!***************************************************************************************************
! This routine takes two lists of integer sequences and compares them to see if they are the
! same. The input lists are assumed to contain unique sequences and perhaps padded with zeros.
FUNCTION lists_match(list1, list2)
integer, intent(in) :: list1(:,:), list2(:,:) 
logical lists_match, rowmatch
integer nL ! length of each list, length of lists
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
! This routine takes a list of label permutations and a list of labelings and reduces the labelings
! list by removing those that are duplicate under the label rotation.
SUBROUTINE remove_label_rotation_dups(lr,lab,labTabin,trgrp,k,d,eps)
integer, intent(in) :: lr(:,:) ! Label rotations (permutations) and temp storage
integer, pointer :: lab(:,:) ! labelings that are unique except for label-rotation permutations
character, intent(in) :: labTabin(:) ! Complete table of labeling "flags", F, N, D, E, etc. 
integer, pointer :: trgrp(:,:) ! translation group, permutations that lead to equivalent labelings
integer, intent(in):: k ! Number of colors in the labelings
integer, intent(in):: d(3) ! Diagonal elements of the SNF
real(dp), intent(in) :: eps ! Finite precision tolerance

character :: labTab(size(labTabin))
integer i, j, n, ip, np, il, nl, idx, ic, nUql, itr, nlr(1), ilr, ia, status
integer b(size(lab,2)), c(size(lab,2)), multiplier(size(lab,2)), ctemp(size(lab,2))
integer trivPerm(size(lab,2)), kc(size(lab,2))
integer, pointer :: perm(:,:) => null()

np = factorial(k)
call get_permutations((/(i,i=0,k-1)/),perm) ! Since we enter this loop so often (for same k) perhaps
! this should be an input

nlr = count(lr/=0,2)
n = size(lab,2)
trivPerm = (/(i,i=1,n)/)
nl = size(lab,1)
multiplier = k**(/(i,i=n-1,0,-1)/)

labTab = labTabin ! Make a copy of the flags table (so that we can cross off duplicates)
! Now that we have a list of label permutations, use them to shrink the input list of labelings
do ilr = 1, nlr(1) ! There's one permutation for each rotation that fixes the lattice
   if (all(lr(:,ilr)==trivPerm)) cycle ! If the permutation is the trivial one, skip over it
   do il = 1, nl ! loop over each labeling in the list
      b = lab(il,:)  ! Get the il-th labeling (we don't want to remove this one, just label-rotated copies of itself)
      do itr = 1,size(trgrp,1) ! Loop over all possible translations of this permutation
         c = b(lr(:,ilr)) ! Permute the labels according to the permutation of this rotation
         c = c(trgrp(itr,:)) ! Permute the labeling according to each translation
         ctemp = c
         do ip = 1,np ! Loop over all permutations (exchanges) of the labels (stored in 'perms')
            do ia=1,k ! Convert the k-ary labeling to one with the labels permuted
               where(ctemp==ia-1); c(:) = perm(ia,ip);endwhere ! b is the permuted labeling
            end do
            if (all(b==c)) cycle ! This permutation didn't change the labeling         
            idx = sum(c*multiplier)+1 ! Find the index for this permutation in the table
            if (labTab(idx)=='F' ) then ! This structure is a label rotation duplicate
               idx = sum(b*multiplier)+1
               labTab(idx) = 'R'
            endif
         enddo
      enddo
   enddo
enddo
! Now replace the input list of labelings with the new list that has shrunk by the number of label-rot dups
! Probably easiest to do this by using the k-nary counter again and just making a list of k-nary n-digit 
! numbers that aren't marked off in the label table.
nUql = count(labTab=='F')
if(associated(lab)) deallocate(lab)
allocate(lab(nUql,n),STAT=status)
if(status/=0) stop "Allocation of lab failed in remove_label_rotation_dups"
kc = 0; ic = 0
do ! Loop over all values of a k-nary, n-digit counter
   idx = sum(kc*multiplier)+1
   if (labTab(idx)=='F') then
      ic = ic + 1   ! Count the number of labelings found so far
      lab(ic,:) = kc ! Store this labeling
   endif
   j = n
   do ! Check to see if we need to roll over any digits, start at the right
      if (kc(j) /= k - 1) exit ! This digit not ready to roll over, exit the loop and advance digit
      kc(j) = 0  ! Rolling over so set to zero
      j = j - 1;       ! Look at the next (going leftways) digit
      if (j < 1) exit  ! If we updated the leftmost digit then we're done
   enddo
   if (j < 1) exit ! We're done counting, exit
   kc(j) = kc(j) + 1 ! Update the next digit (add one to it)
enddo
if (ic/=nUql) stop "relabeling error"

ENDSUBROUTINE remove_label_rotation_dups

!****************************************************************************************************
! This routine takes in the size of the three cyclic groups (i.e., the diagonal elements of the
! Smith normal form (SNF) and generates all labelings that are unique. That is, it eliminates
! "translation" duplicates as well as superperiodic duplicates (non-primitive
! superstructures). It now also removes "label-permutation" duplicates---labelings that are not
! unique when the labels themselves (not their positions) are permuted (e.g., 00111 <-->
! 11000).  The basic idea of the routine is to run like an "odometer", generating all numbers
! (base k) from 0 to k^n - 1. Then use the translations (permutations of the group) to
! eliminate translation duplicates and superperiodic cases.  ** We must also eliminate cases
! where the superlattice is fixed by a parent lattice symmetry but the labels are
! permuted. These cases are also duplicates.  But this is best done outside of and after this
! subroutine because those duplicates are unique to the HNF. **

SUBROUTINE generate_labelings(k,d,l,lab,trgrp,full)
integer, intent(in) :: k,d(3) ! Number of colors/labels, SNF diagonal elements
integer, pointer :: l(:,:) ! labelings
character, pointer :: lab(:) ! Array to store markers for every raw labeling
! E=>missing digits, F=>unique, D=>translation duplicate, N=>non-primitive, R=>label-rotation dup
! Need to pass lab out so that it can be used to remove label-rotation duplicates later
integer, pointer :: trgrp(:,:) ! array for storing the translation group. Need to pass this out
! later too for when the label-rotation duplicates are removed.
logical, intent(in) :: full ! specify whether the full labelings list should be used or not

integer cnt ! Number of unique labelings (double check on the generator)
integer j ! Index variable (place index) for the k-ary counter
integer ic, i, q, ia ! loop counters, index variables
integer nexp ! number of raw labelings that the k-ary counter should generate
integer n ! number of labels in each labeling (i.e., determinant size)
integer idx ! the base 10 equivalent of the current base k labeling
integer a(d(1)*d(2)*d(3)), b(d(1)*d(2)*d(3)) ! the "odometer"; label-permuted odometer
integer multiplier(d(1)*d(2)*d(3)) ! place values for each digit k^(i-1) for the i-th digit
integer c(0:k-1) ! running sum (count) of the number of each label type 
integer id, iq ! Counter for labels that are duplicates, for those unique
integer, allocatable :: tl(:,:) ! temporary storage for output variable l (labelings)
integer, pointer :: perms(:,:) ! List of permutations of the k labels
integer :: np, ip, status ! Loops over permutations, allocate error flag

if (associated(lab)) deallocate(lab)
allocate(lab(k**(d(1)*d(2)*d(3))),STAT=status)
if(status/=0) stop "Failed to allocate memory for 'lab' in generate_labelings"

n = product(d); nexp = k**n  ! Number of digits in k-ary counter; upper limit of k-ary counter
a = 0; multiplier = k**(/(i,i=n-1,0,-1)/) ! The counter; multiplier to convert to base 10
lab = ''; iq = 0  ! Index for labelings; number of unique labelings
if (k>12) stop "Too many labels in 'generate_labelings'"

np = factorial(k) ! Number of permutations of labels (not labelings)
call get_permutations((/(i,i=0,k-1)/),perms)
call count_full_colorings(k,d,cnt,full) ! Use the Polya polynomial to count the labelings
call make_translation_group(d,trgrp) ! Find equivalent translations (permutations of labelings)
allocate(tl(cnt,n),STAT=status)   ! Temp storage for labelings
if(status/=0) stop "Failed to allocate memory for tl in generate_labelings"

ic = 0; c = 0; c(0) = n ! Loop counter for fail safe; initialize digit counter
do; ic = ic + 1
   if (ic > nexp) exit ! Fail safe
   idx = sum(a*multiplier)+1;  ! Index of labeling in base 10
   if (idx/=ic) stop "index bug!"
   if (any(c==0)) then ! Check to see if there are missing digits
      id = id + 1; ! Keep track of the number of incomplete labelings
      if (lab(idx)=='' .and. .not. full) then ! If it isn't marked and we want a partial list, mark it as "error"
         lab(idx) = 'E';    ! Could mark its brothers too...
      endif
   endif
   ! If this label hasn't been marked yet, mark it as found 'F'
   ! and mark its duplicates as 'D'
   if (lab(idx)=='') then
      lab(idx) = 'F'
      do q = 2,n ! Mark translation duplicates and eliminate superperiodic (non-primitive) colorings
         idx = sum(a(trgrp(q,:))*multiplier)+1
         if (idx==ic) lab(idx)='N' ! This will happen if the coloring is superperiodic
                                   ! (i.e., non-primitive superstructure)
         if (lab(idx)=='') lab(idx) = 'D'  ! Mark as a translation duplicate
      enddo
      do q = 1,n ! Loop over translates; Mark label-permutation duplicates ("symmetric" concentrations)  
         do ip = 2,np ! Loop over all permutations of the labels (stored in 'perms')
            forall(ia=1:k) ! Convert the k-nary labeling to one with the labels permuted
               where(a==ia-1); b(:) = perms(ia,ip);endwhere ! b is the permuted labeling
            end forall
            idx = sum(b(trgrp(q,:))*multiplier)+1
            if  (lab(idx)=='') lab(idx) = 'N'
         enddo
      enddo
      ! Store the unique labelings
      if (lab(ic)=='F') then ! this one isn't superperiodic (it is primitive) and it is unique
         iq = iq + 1
         tl(iq,:) = a
      endif
   endif
   ! Advance the counter and keep track of the # of numbers
   j = n ! Reset the digit index (start all the way to the right again)
   do ! Check to see if we need to roll over any digits, start at the right
      if (a(j) /= k - 1) exit ! This digit not ready to roll over, exit the loop and advance digit
      a(j) = 0  ! Rolling over so set to zero
      c(k-1) = c(k-1) - 1  ! Update the digit count; just lost a digit of type k-1 (largest possible)
      c(0) = c(0) + 1  ! So we pick up another zero in the current place ([k-1]->0)
      j = j - 1;       ! Look at the next (going leftways) digit
      if (j < 1) exit  ! If we updated the leftmost digit then we're done
   enddo
   if (j < 1) exit ! We're done counting, exit
   a(j) = a(j) + 1 ! Update the next digit (add one to it)
   c(a(j)) = c(a(j)) + 1 ! Add 1 to the number of digits of the j+1-th kind
   c(a(j)-1) = c(a(j)-1) - 1     ! subtract 1 from the number of digits of the j-th kind
   if (sum(c) /= n .and. .not. full) stop 'counting bug'
enddo
if (ic /= nexp) stop 'Bug: Found the wrong number of labels!'
! Store the results
if(associated(l)) deallocate(l)
allocate(l(iq,n),STAT=status)
if(status/=0) stop "Allocation of 'l' failed in generate_labelings"
l = tl(1:iq,:)
END SUBROUTINE generate_labelings
 
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
! to get another member of the group. This list of mappings in the list of labelings (colorings) 
! that leave the superstructure unchanged.
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
   do i = 1,n ! This approach is an N^2 loop. Can this be improved? Does it matter?
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
