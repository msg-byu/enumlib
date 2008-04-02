MODULE combinatorics
use num_types
use crystal_types
use vector_matrix_utilities
use numerical_utilities
use rational_mathematics, only: gcd
implicit none
private
public  remove_label_rotation_dups, get_permutations, m_partitions_of_n, factorial, count_full_colorings, &
       generate_labelings, nchoosek, k_ary_counter, make_member_list
INTERFACE factorial
   MODULE PROCEDURE factorial_int, factorial_rank1
END INTERFACE
CONTAINS
!***************************************************************************************************
! This routine takes a list of SNF left transformation matrices and the parent lattice and finds
! a list of permutations (label-rotation permutations, which leave the superlattice itself fixed).
! The routine then takes a list of labelings and removes the duplicates
SUBROUTINE remove_label_rotation_dups(L,A,R,G,lab,labTabin,trgrp,k,d,eps)
integer, intent(in) :: L(:,:) !Left transformation matrices for SNF (input)
integer, intent(in) :: G(:,:) ! image group (input)
integer, pointer :: lab(:,:) ! labelings that are unique except for label-rotation permutations
real(dp) :: A(3,3) ! The parent lattice
!type(opList), pointer :: R(:)
real(dp) :: R(:,:,:)
character, intent(in) :: labTabin(:)
integer, pointer :: trgrp(:,:)
integer k ! Number of colors in the labelings
integer d(3) ! Diagonal elements of the SNF

integer, allocatable :: tl(:,:) ! temporary list of labels
integer, allocatable :: Gp(:,:)  ! G prime, the permuted (by T) image group
character :: labTab(size(labTabin))
integer i, j, n, iRot, nRot, il, nl, idx, ic, nUql, itr, np, ip, ia, allocStat
logical err
real(dp), dimension(3,3) :: T, Ident, A1, A1inv, invA
integer lr(size(lab,2),48)
integer b(size(lab,2)), c(size(lab,2)), multiplier(size(lab,2)), ctemp(size(lab,2))
real(dp) :: eps
integer trivPerm(size(lab,2)), kc(size(lab,2))
integer, pointer :: perm(:,:)

np = factorial(k)
call get_permutations((/(i,i=0,k-1)/),perm)
!print *, labtabin
!write(*,'("<<< Count",i3)') count(labtabin=='F')
!do i = 1,3
!   do j =1,size(R,3)
!      write(*,'(3f6.2,3x)',advance = 'no') R(i,:,j)
!      enddo; print *
!enddo
!print *
n = size(lab,2)
trivPerm = (/(i,i=1,n)/)
nl = size(lab,1)
multiplier = k**(/(i,i=n-1,0,-1)/)
Ident = 0; Ident(1,1) = 1;Ident(2,2) = 1;Ident(3,3) = 1 ! Identity matrix
!nHNF = size(L,3) ! The number of HNF matrices
allocate(tl(size(lab,1),size(lab,2)),STAT=allocStat) ! Allocate the temporary list
if(allocStat/=0) stop "Allocation problem in remove_label_ratation_dups: tl"
allocate(Gp(size(G,1),size(G,2)))
labTab = labTabin

!do i = 1,size(lab,1)
!   write(*,'(40i1)') lab(i,:)
!enddo

tl = lab ! Copy the labels
call matrix_inverse(A,invA,err)  ! Need A^-1 to form the transformation

! Find the list of permutations that correspond to label-rotations
A1 = matmul(L,invA)
call matrix_inverse(A1,A1inv,err)
nRot = size(R,3)
lr = 0 ! Reset the list of label-rotation permutations
do iRot = 1, nRot ! Use each rotation that fixes the lattice.
!   write(*,'("Rot:",i3)') iRot
   if (equal(R(:,:,iRot),Ident,eps)) cycle ! Skip the identity---it won't rotate the labels
!   do i = 1,3
!      write(*,'("Rot: ",3f7.3)') R(i,:,iRot)
!   enddo;print *
   T = matmul(A1, matmul(R(:,:,iRot),A1inv))  ! This is the transformation
   if (.not. equal(T,nint(T),eps)) then
      print *, 'ERROR: remove_label_rotation_dups: Transformation is not integer'
      stop
   endif
!   do i = 1,3
!      write(*,'("T : ",20i3)') nint(T(i,:))
!   enddo;print *
   ! Now need to use nint(T) to rearrange the image group
   Gp = matmul(nint(T),G)  ! Gp (prime) is the tranformed image group. Need this to get the permutation
!   do i = 1,3
!      write(*,'("Gp: ",20i3)') Gp(i,:)
!   enddo; print *

   do i=1,3; Gp(i,:) = modulo(Gp(i,:),d(i));enddo
!      do i = 1,3
!         write(*,'("G : ",20i3)') G(i,:)
!      enddo;print *
!      do i = 1,3
!         write(*,'("Gp: ",20i3)') Gp(i,:)
!      enddo; print *
!
   do i = 1, n
      do j = 1, n
         if (all(Gp(:,j)==G(:,i))) then ! the two images are the same for the i-th member
            lr(i,iRot) = j
         endif
      enddo
   enddo
   if (any(lr(:,iRot)==0)) stop "Transform didn't work. Gp is not a permutation of G"
enddo

! Now that we have a list of label permutations, use them to shrink the input list of labelings
do iRot = 1, nRot ! There's one permutation for each rotation that fixes the lattice
!   print *
!   write(*,'("Rotation perm: ",i3)') iRot
!   write(*,'("perm: ",40i2)') lr(:,iRot)
   if (all(lr(:,iRot)==trivPerm)) cycle ! If the permutation is the trivial one, skip over it
   if (all(lr(:,iRot)==0)) cycle ! This happens for the trivial rotation (identity)
   do il = 1, nl ! loop over each labeling in the list
      b = lab(il,:)  ! Get the il-th labeling (we don't want to remove this one,
      ! just label-rotated copies of itself
!      print *,"labeling to rotate:",b
      do itr = 1,size(trgrp,1) ! Loop over all possible translations of this permutation
         c = b(lr(:,iRot)) ! Permute the labels according to the permutation of this rotation
         c = c(trgrp(itr,:)) ! Permute the labeling according to each translation
         ctemp = c
         do ip = 1,np ! Loop over all permutations of the labels (stored in 'perms')
            do ia=1,k ! Convert the k-ary labeling to one with the labels permuted
               where(ctemp==ia-1); c(:) = perm(ia,ip);endwhere ! b is the permuted labeling
            end do
            if (all(b==c)) cycle ! This permutation didn't change the labeling         
            idx = sum(c*multiplier)+1 ! Find the index for this permutation in the table
  !          write(*,'(a1,":",i4,4x,40i2)') labtab(sum(b*multiplier)+1),idx,b
   !         write(*,'(a1,":",i4,4x,40i2)') labtab(idx),idx,c
            if (labTab(idx)=='F' ) then
               idx = sum(b*multiplier)+1
 !              print *, "found a rotation duplicate"
               labTab(idx) = 'R'
            endif
         enddo
      enddo
   enddo
enddo
! Now replace the input list of labelings with the new list that has shrunk by the number of label-rot dups
! Probably easiest to do this by using the k-ary counter again and just making a list of k-ary n-digit 
! numbers that aren't marked off in the label table.
nUql = count(labTab=='F')
!print *, labTab
!write(*,'("<<< after count",i3)') nUql
if(associated(lab)) deallocate(lab)
allocate(lab(nUql,n),STAT=allocStat)
if(allocStat/=0) stop "Allocation problem in remove_label_ratation_dups: lab")
kc = 0; ic = 0
do ! Loop over all values of a k-ary, n-digit counter
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
!do i = 1,nUql
!   write(*,'(40i2)') lab(i,:)
!enddo


ENDSUBROUTINE remove_label_rotation_dups

!****************************************************************************************************
! This routine takes in the size of the three cyclic groups (i.e., the diagonal elements of the
! Smith normal form (SNF) and generates all labelings that are unique. That is, it eliminates
! "translation" duplicates as well as superperiodic duplicates (non-primitive superstructures). It now
! also removes "label-permutation" duplicates---labelings that are not unique
! when the labels themselves (not their positions) are permuted (e.g., 00111 <--> 11000).
!   The basic idea of the routine is to run like an "odometer", generating all numbers (base k) from
! 0 to k^n - 1. Then use the translations (permutations of the group) to eliminate translation
! duplicates and superperiodic cases.  ** We must also eliminate cases where the superlattice is fixed
! by a parent lattice symmetry but the labels are permuted. These cases are also duplicates.  But this is 
! best done outside of and after this subroutine. **
SUBROUTINE generate_labelings(k,d,l,lab,trgrp)
integer, intent(in) :: k,d(3) ! Number of colors/labels, SNF diagonal elements
integer, pointer :: l(:,:) ! labelings
character, pointer :: lab(:) ! Array to store markers for every raw labeling
! E=>missing digits, F=>unique, D=>translation duplicate, N=>non-primitive, R=>label-rotation dup
! Need to pass lab out so that it can be used to remove label-rotation duplicates later
integer, pointer :: trgrp(:,:) ! array for storing the translation group. Need to pass this out
! later too for when the label-rotation duplicates are removed.

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
integer :: np, ip ! Loops over permutations
integer allocStat
if (associated(lab)) deallocate(lab)
allocate(lab(k**(d(1)*d(2)*d(3))),STAT=allocStat) 
if(allocStat/=0) stop "Allocation problem in generate_labelings: lab" 

n = product(d); nexp = k**n  ! Number of digits in k-ary counter; upper limit of k-ary counter
a = 0; multiplier = k**(/(i,i=n-1,0,-1)/) ! The counter; multiplier to convert to base 10
lab = ''; iq = 0  ! Index for labelings; number of unique labelings
if (k>12) stop "Too many labels in 'generate_labelings'"

np = factorial(k) ! Number of permutations of labels (not labelings)
call get_permutations((/(i,i=0,k-1)/),perms)
call count_full_colorings(k,d,cnt) ! Use the Polya polynomial to count the labelings
call make_translation_group(d,trgrp) ! Find equivalent translations (permutations of labelings)


! Temp storage for labelings
allocate(tl(cnt,n),STAT=allocStat) ! Allocate the temporary list
if(allocStat/=0) stop "Allocation problem in generate_labelings: tl"  

ic = 0; c = 0; c(0) = n ! Loop counter for fail safe; initialize digit counter
do; ic = ic + 1
   if (ic > nexp) exit ! Fail safe
   idx = sum(a*multiplier)+1;  ! Index of labeling in base 10
   if (idx/=ic) stop "index bug!"
   if (any(c==0)) then ! Check to see if there are missing digits
      id = id + 1;
      if (lab(idx)=='') then ! If it isn't marked, mark it as "error"
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
!         write(*,'("trans: ",9i1,4x,a1)') a(trgrp(q,:)),lab(idx) 
         !write(*,*) "trans:", a(trgrp(q,:)),lab(idx) 
      enddo
      do q = 1,n ! Loop over translates; Mark label-permutation duplicates ("symmetric" concentrations)  
         do ip = 2,np ! Loop over all permutations of the labels (stored in 'perms')
            forall(ia=1:k) ! Convert the k-ary labeling to one with the labels permuted
               where(a==ia-1); b(:) = perms(ia,ip);endwhere ! b is the permuted labeling
            end forall
            idx = sum(b(trgrp(q,:))*multiplier)+1
            if  (lab(idx)=='') lab(idx) = 'N'
 !           write(*,'("labpr: ",9i1,4x,a1)') b(trgrp(q,:)),lab(idx)
            !write(*,*) "labpr:",b(trgrp(q,:)),lab(idx)
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
   if (sum(c) /= n) stop 'counting bug'
enddo
if (ic /= nexp) stop 'Bug: Found the wrong number of labels!'
! Store the results
if(associated(l)) deallocate(l)
allocate(l(iq,n)STAT=allocStat) 
if(allocStat/=0) stop "Allocation problem in generate_labelings: l"
l = tl(1:iq,:)
END SUBROUTINE generate_labelings
 
!**************************************************************************************************
! Takes the length of three cyclic group and constructs the member list so that
! the permutations can be determined. Each member has three components, corresponding to the
! entries for each of the three cyclic groups. 
SUBROUTINE make_member_list(n,p)
integer, intent(in)  :: n(3)! Diagonal elements of the SNF
integer, pointer :: p(:,:)  ! List of members of the translation group

integer im  ! loop over members
if (associated(p)) deallocate(p)
allocate(p(3,product(n)))
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
!do im = 1,3;write(*,'(20i2)') p(im,:);enddo
ENDSUBROUTINE make_member_list

!***************************************************************************************************
! This routine finds all the permutations of the group members that leave the decoration unchanged.
! Essentially we are finding a list of mappings: add to the group one of the members of the group
! to get another member of the group. This list of mappings in the list of labelings (colorings) 
! that leave the superstructure unchanged.
SUBROUTINE make_translation_group(d,trans)
integer, intent(in) :: d(3) ! members of the group, diagonal elements of SNF
integer, pointer :: trans(:,:) ! Translations that leave the superstructure unchanged

integer n, im,i, j
integer tg(3,d(1)*d(2)*d(3)) ! Temporary storage for the translations
integer, pointer :: m(:,:) ! List of group members

! Get a list of group members; essentially a mixed-radix, 3 digit counter
call make_member_list(d,m)
n = product(d) ! Number of elements in each permutation group
if (associated(trans)) deallocate(trans)
allocate(trans(n,n))

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
SUBROUTINE count_full_colorings(k,d,count)
integer, intent(in) :: k, d(3)  ! Number of colors, diagonal entries of the Smith Normal Form matrix
integer, intent(out) :: count   ! Number of unique, full colorings

integer x1,x2,x3 ! Counters over d1,d2,d3, in the triple sum
integer p        ! Counter over number of terms in the exclusion/inclusion counting
integer m, tc    ! Color index (counter), intermediate (temporary) counter
integer n        ! Number of elements to be colored
m = k; count = 0
n = product(d)
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
ENDSUBROUTINE count_full_colorings

!***************************************************************************************************
FUNCTION nchoosek(n,k)
integer nchoosek
integer, intent(in) :: n,k
integer i
if (k > n) stop "Error in arguments to nchoosek"
if (k > n/2) then
   nchoosek = product( (/ (i,i=n,k+1,-1) /) )/factorial(n-k)
else
   nchoosek = product( (/ (i,i=n,n-k+1,-1) /) )/factorial(k); endif
ENDFUNCTION nchoosek

!*******************************************************************************
! This subroutine generates all permutations of the list. Swapping identical
! items does *not* give another permutation. List is generated in
! lexicographical order. Input list must be in order (or the routine quits with 
! an error)
!
! This could be improved by determining a priori the number of matching elements
! and thus the number of actual permuations. Otherwise, the temp storage might
! get too big too fast.
SUBROUTINE get_permutations(list,perms)
integer, intent(in) :: list(:)
integer, pointer :: perms(:,:)

integer, allocatable:: tempPerms(:,:)
integer Nitems, Np, a(size(list)),i
integer j, l, temp, ip

! If Nitems > 12 or so, this is a bad way of doing things...need to find a stable way of calculating a binomial, trinomial, etc. 

Nitems = size(list)
Np = factorial(Nitems)  ! This is the upper limit for the number of permutations
Np = 1
do i = Nitems, Nitems/2, -1
   Np = Np * i
enddo
!print *, Np, Np/2
Np = Np/factorial(Nitems/2)
!write(*,'("permutations. Nitems: ",i6,"Num. perms: ",i16)') Nitems, Np

!if(allocated(tempPerms)) deallocate(tempPerms)
allocate(tempPerms(Nitems,Np),STAT=allocStat) 
if(allocStat/=0) stop "Allocation problem in get_permutations: tempPerms"

a = list ! Copy of the original input list
! Make sure the items are initially ordered
do i = 1, Nitems-1
   if (a(i+1)<a(i)) stop "ERROR in 'get_permutations'---incorrectly ordered input"
enddo
ip = 0
outer:do
    ip = ip + 1;
    tempPerms(:,ip)=a
    do j = Nitems-1,0,-1
        if (j==0) exit outer
        if (a(j) < a(j+1)) exit
    enddo
    do l = Nitems,1,-1
        if (a(j) < a(l)) exit
    enddo
    temp = a(j)
    a(j) = a(l)
    a(l) = temp
    a(j+1:Nitems)=a(Nitems:j+1:-1)
enddo outer
if(associated(perms)) deallocate(perms)     
allocate(perms(Nitems,ip),STAT=allocStat) 
if(allocStat/=0) stop "Allocation problem in get_permutations: perms"

perms = tempPerms(:,:ip) ! There will be fewer permutations than N! if there
deallocate(tempPerms)    ! are identical items in the permuted list
END SUBROUTINE get_permutations

!***************************************************************************************************
! This routine generates a list of all numbers between 0 and n^k-1
SUBROUTINE k_nary_counter(list,base,n)
integer, pointer :: list(:,:)
integer, intent(in) :: base 
integer, intent(in) :: n ! number of digits
integer k, j, il, allocStat
integer :: a(n)

k = base; allocate(list(n,k**n),STAT=allocStat) 
if(allocStat/=0) stop "Allocation problem in k_nary_counter: list"

il = 1; a = 0
list(:,il) = a
do
   j = n ! Number of digits
   do ! Check to see if we need to roll over any digits, start at the right
      if (a(j) /= k - 1) exit 
      a(j) = 0  ! Rolling over so set to zero
      j = j - 1;       ! Look at the next (going leftways) digit
      if (j < 1) exit  ! If we updated the leftmost digit then we're done
   enddo
   if (j < 1) exit ! We're done counting, exit
   a(j) = a(j) + 1 ! Update the next digit
   il = il + 1
   list(:,il) = a
enddo
ENDSUBROUTINE k_nary_counter

!*******************************************************************************
! This subroutine generates all partions of n into m blocks. This is needed to
! generate all possible concentration vectors for each structures (then within
! each structure we need all permutations of atomic configurations...)
! See page 38 in Knuth V. 4 Fascicle 3
SUBROUTINE m_partitions_of_n(n,m,part)
integer, intent(in) :: n,m ! number to be partitioned, number of block
integer, pointer   :: part(:,:)

integer :: a(m)
integer          :: ip  ! counter for number of partitions
integer          :: x, s ! Temporary variables for swapping/reassigning values
integer          :: j    ! index pointer for next position to be incremented

! First we need to count the number of partitions to know how big to allocate
! the "part" array (is there a closed form expression for this? I don't think so
! but...) 
if (m>n) then; stop "Bad input for 'm_partitions_of_n' routine";endif
if (m<2) then; stop "Trivial case 'm_partitions_of_n' routine";endif
! Initialize the first partition
a(1) = n - (m-1); a(2:m) = 1
ip = 0
do; ip = ip + 1;
    if (a(2) < a(1)-1) then ! "crumble" at left if possible
        a(1) = a(1)-1
        a(2) = a(2)+1
        cycle; endif
    if (m==2) exit ! For binary case, we're done here
    ! Find the leftmost position that can be increased
    j = 3 ! This is the leftmost possibility but...
    s = a(1)+a(2)-1
    do while (a(j)>=a(1)-1)
        s = s+a(j)
        j = j+1
        if (j>m) exit
    enddo ! Now s = part(1)+part(2)+...+part(j-1)-1)
    if (j > m) exit ! We're done counting
    x = a(j)+1 ! Store the value of the j-th position incremented by one
    a(j) = x   ! Make the j-th position have this value
    j = j-1    ! Now look one to the left
    do while (j>1)
        a(j) = x  ! Make this next-left position match the slot that was updated
        s = s - x ! Keep track of how many of the original 'n' are left
        j = j-1   ! Now look to the left again...
    enddo
    a(1) = s  ! Now make the leftmost slot take all the leftover...
enddo
!if (allocated(part)) deallocate(part)
allocate(part(ip,m),STAT=allocStat) 
if(allocStat/=0) stop "Allocation problem in m_partitions_of_n: part"

ip = 0; a(1) = n - (m-1); a(2:m) = 1
do; ip = ip + 1
    ! Visit the partition
    part(ip,:) = a ! store the ip-th partition
    if (a(2) < a(1)-1) then ! "crumble" at left if possible
        a(1) = a(1)-1
        a(2) = a(2)+1
        cycle; endif
    if (m==2) exit ! For binary case, we're done here
    ! Find the leftmost position that can be increased
    j = 3 ! This is the leftmost possibility but...
    s = a(1)+a(2)-1
    do while (a(j)>=a(1)-1)
        s = s+a(j)
        j = j+1
        if (j>m) exit
    enddo ! Now s = part(1)+part(2)+...+part(j-1)-1)
    if (j > m) exit ! We're done counting
    x = a(j)+1
    a(j) = x
    j = j-1
    do while (j>1)
        a(j) = x
        s = s - x
        j = j-1
    enddo
    a(1) = s
enddo
END SUBROUTINE m_partitions_of_n

FUNCTION factorial_int(N)
integer N,i,factorial_int
factorial_int = 1
do i=2,N; factorial_int=factorial_int*i;enddo
END FUNCTION factorial_int

FUNCTION factorial_rank1(N)
integer N(:),i,factorial_rank1(size(N)),j
factorial_rank1 = 1
do j=1,size(N)
   do i=2,N(j); factorial_rank1(j)=factorial_rank1(j)*i;enddo;enddo
END FUNCTION factorial_rank1

END MODULE combinatorics
