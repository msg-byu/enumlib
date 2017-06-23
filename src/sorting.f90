MODULE sorting
use num_types

implicit none
private
public heapsort, sort_permutations_list

INTERFACE heapsort
   MODULE PROCEDURE heapsort_records_int, heapsort_records_dp, sort_concs
END INTERFACE
CONTAINS

  !!<summary>This routine takes a list of integers and sorts them from
  !!least to greatets. It also sorts a second list of labels so that the
  !!relative order of the entries of the two lists remains the same.</summary>
  !!<parameter name="list" regular="true">The 1D integer array to be
  !!sorted.</parameter>
  !!<parameter name="labels" regular="true">The list of labels for the
  !!concentrations.</parameter>
  !!<parameter name="conc_map" regular="true">The mapping of arrowed
  !!to non arrowed concentrations</parameter>
  SUBROUTINE sort_concs(list,labels,conc_map)
    integer, intent(inout) :: list(:), labels(:)
    integer, intent(in) :: conc_map(:,:)
    !!<local name="sorted_list">A temporary array for storing the sorted array.</local>
    !!<local name="minl">The location of the smallest entry in the list.</local>
    !!<local name="species_i">Counting index.</local>
    !!<local name="species_j">Counting index.</local>
    !!<local name="species_a">Counting index.</local>
    !!<local name="arrow_labels">Temporary array of the arrow_labels.</local>
    !!<local name="arrow_concs">Temporary arroy of the arrow_concs</local>
    integer :: sorted_list(size(list)), sorted_labels(size(labels)), arrow_labels(size(labels)),arrow_concs(size(list))
    integer :: minl
    integer :: species_i, species_j, species_a

    ! If there are no arrows present then we can just sort the
    ! lists like normal.
    if (all(conc_map == 0)) then
       do species_i =1, size(list)
          minl = minloc(list, 1, list .NE. -1)
          sorted_list(species_i) = list(minl)
          sorted_labels(species_i) = labels(minl)
          list(minl) = -1
       end do
    else
       ! If there are arrows present then we want to make sure that
       ! they appear at the end of the concentration list. We do this
       ! by sorting the arrows colors and the non arrow colors into
       ! two seperate lists then joining them at the end.
       arrow_concs = -1
       arrow_labels = 0
       species_j = 1
       species_a = 1
       do species_i =1, size(list)
          minl = minloc(list, 1, list /= -1)
          if (any(conc_map(:,1) == labels(minl)+1)) then
             ! This is an arrowed species, store it in the arrowed list.
             arrow_concs(species_a) = list(minl)
             arrow_labels(species_a) = labels(minl)
             list(minl) = -1
             species_a = species_a + 1
          else
             ! store the species in the unarrowed list.
             sorted_list(species_j) = list(minl)
             sorted_labels(species_j) = labels(minl)
             list(minl) = -1
             species_j = species_j + 1
          end if
       end do
       ! Now we place the arrows at the end of the sorted list.
       do species_i =1, species_a -1
          sorted_list(species_j) = arrow_concs(species_i)
          sorted_labels(species_j) = arrow_labels(species_i)
          species_j = species_j + 1
       end do
    end if
    labels = sorted_labels
    list = sorted_list

  END SUBROUTINE sort_concs
    

  !!<summary>This routine takes a list of integer sequences (rows) and
  !!removes duplicates so that each row is unique. It seems that the
  !!check for uniqueness can be done in O(N) time if the list is in
  !!"alphabetical" order. So it's probably most efficient, at least
  !!for big lists, to sort it first (using an efficent sort), and then
  !!make the unique list.</summary>
  !!<parameter name="key">The values to sort</parameter>
  SUBROUTINE sort_permutations_list(key)
    !integer, pointer :: list(:)! indices of the records to be
    ! rearranged using the keys (R in Knuth) Not going to use the
    ! list, just sort the keys
    integer, pointer:: key(:,:)

    integer N, l, r, i, j ! Same as in Knuth
    integer pK(size(key,2)) ! "plain" K (c'mon Knuth, bad notation!)
    integer Ng, tkey(size(key,1),size(key,2)), count ! Sequence
    ! length, temp. key, # of unique sequences
    logical duplicate

    N = size(key,1)
    Ng = size(key,2) ! Number of elements in each sequence (number of group members)

    l = N/2 + 1; r = N; 
    do ! H2: decrease l or r
       if (l>1) then 
          l = l - 1; pK = key(l,:)
       else
          pK = key(r,:); key(r,:) = key(1,:)
          r = r - 1
          if (r==1) then; key(1,:) = pK; exit; endif
       endif
       ! H3: prepare for siftup
       j = l 
       ! H4: Advance downward
       do; i = j; j = 2*j
          if (j<r) then ! H5: Find larger child
             if (greater_than(key(j+1,:),key(j,:))) j = j + 1
             if (greater_than(pK,key(j,:))) exit ! Go to step H8
          else if(j==r) then! H6: Larger than pK?
             if (greater_than(pK,key(j,:))) exit ! Go to step H8
          else
             exit ! Go to step H8
          endif
          ! H7: Move it up
          key(i,:) = key(j,:)
       enddo
       ! H8: Store key
       key(i,:) = pK
    enddo ! Return to step H2
    do i = 1, N-1  ! Double check that sorting really worked...
       if (greater_than(key(i,:), key(i+1,:))) stop "Sorting routine (sort_permutations_list) failed"
    enddo

    ! Now reduce the list "key" to a list without duplicates. This is
    ! just an O(N) operation since the list is sorted (i.e., the
    ! duplicates must be neighbors in the list)
    count = 1
    tkey(1,:) = key(1,:)
    do i = 2, N
       duplicate = .false.
       if (all(key(i,:)==key(i-1,:))) then
          duplicate = .true.; cycle
       endif
       if (.not. duplicate) then
          count = count + 1
          tkey(count,:) = key(i,:)
       endif
    enddo
    deallocate(key)
    allocate(key(count,Ng))
    key = tkey(1:count,:)
    
  CONTAINS
    !!<summary></summary>
    !!<parameter name="vecA" regular="true"></parameter>
    !!<parameter name="vecB" regular="true"></parameter>
    FUNCTION greater_than(vecA, vecB)
      integer, intent(in):: vecA(:), vecB(:)
      logical :: greater_than
      integer :: i

      greater_than = .false.
      do i = 1, size(vecA)
         if(vecA(i) > vecB(i)) then; greater_than = .true.
         elseif(vecB(i) > vecA(i))then; exit ! Early exit if the two vectors aren't equivalent
         endif
      enddo
    ENDFUNCTION greater_than
  ENDSUBROUTINE sort_permutations_list


  !!<summary>Uses heapsort to order a group of "records" (list) using
  !!a set of "keys" (the values to use in sorting). See Knuth "Sorting
  !!and Searching" pgs. 142--147. In Knuth's discussion, rearrangments
  !!of the list also imply rearranging the keys. This is not explicit
  !!in Knuth's discussion.</summary>
  !!<parameter name="list">indices of the records to be rearranged
  !!using the keys (R in Knuth)</parameter>
  SUBROUTINE heapsort_records_int(list,key)
    integer, pointer :: list(:)
    integer, intent(inout):: key(:)! values on which to sort, the "keys" (K in Knuth)

    integer :: N, l, r, i, j ! Same as in Knuth
    integer :: pK, pR ! "plain" K, R (c'mon Knuth, bad notation!)

    N = size(key)
    if (associated(list)) deallocate(list)
    allocate(list(N))
    list = (/(i,i=1,N)/) ! Start with an unsorted set of indices
    l = N/2 + 1; r = N; list = (/(i,i=1,N)/)
    do ! H2: decrease l or r
       if (l>1) then 
          l = l - 1; pR = list(l); pK = key(l)
       else
          pR = list(r); pK = key(r); list(r) = list(1); key(r) = key(1)
          r = r - 1
          if (r==1) then; list(1) = pR; key(1) = pK; exit; endif
       endif
       ! H3: prepare for siftup
       j = l 
       ! H4: Advance downward
       do; i = j; j = 2*j
          if (j<r) then ! H5: Find larger child
             if (key(j) < key(j+1)) j = j + 1
             if (pK > key(j)) exit ! Go to step H8
          else if(j==r) then! H6: Larger than pK?
             if (pK > key(j)) exit ! Go to step H8
          else
             exit ! Go to step H8
          endif
          ! H7: Move it up
          list(i) = list(j)
          key(i) = key(j)
       enddo
       ! H8: Store R
       list(i) = pR
       key(i) = pK
    enddo ! Return to step H2
    do i = 1, N-1  ! Double check that sorting really worked...
       if (key(i) > key(i+1)) stop "Sorting routine (heapsort_records_int) failed"
    enddo

  ENDSUBROUTINE heapsort_records_int

  !!<summary>See above except keys are dp reals instead of integers
  !!and see Knuth "Sorting and Searching" pgs. 142--147.</summary>
  !!<parameter name="list">indices of the records to be rearranged
  !!using the keys (R in Knuth)</parameter>
  !!<parameter name="key" regular="true">values on which to sort, the
  !!"keys" (K in Knuth)</parameter>
  SUBROUTINE heapsort_records_dp(list,key) 
    integer, pointer :: list(:)
    real(dp), intent(inout):: key(:)

    integer N, l, r, i, j ! Same as in Knuth
    integer pR ! "plain" R (c'mon Knuth, bad notation!)
    real(dp) pK! "plain" K
    
    N = size(key)
    if (associated(list)) deallocate(list)
    allocate(list(N))
    l = N/2 + 1; r = N; list = (/(i,i=1,N)/)
    do ! H2: decrease l or r
       if (l>1) then 
          l = l - 1; pR = list(l); pK = key(l)
       else
          pR = list(r); pK = key(r); list(r) = list(1); key(r) = key(1)
          r = r - 1
          if (r==1) then; list(1) = pR; key(1) = pK; exit; endif
       endif
       ! H3: prepare for siftup
       j = l 
       ! H4: Advance downward
       do; i = j; j = 2*j
          if (j<r) then ! H5: Find larger child
             if (key(j) < key(j+1)) j = j + 1  ! This might not be a
             ! "stable sort" if we don't account for finite precision
             ! here...
             if (pK > key(j)) exit ! Go to step H8
          else if(j==r) then! H6: Larger than pK?
             if (pK > key(j)) exit ! Go to step H8
          else
             exit ! Go to step H8
          endif
          ! H7: Move it up
          list(i) = list(j)
          key(i) = key(j)
       enddo
       ! H8: Store R
       list(i) = pR
       key(i) = pK
    enddo ! Return to step H2
    do i = 1, N-1  ! Double check that sorting really worked...
       if (key(i) > key(i+1)) stop "Sorting routine (heapsort_records_dp) failed"
    enddo
  ENDSUBROUTINE heapsort_records_dp

END MODULE
