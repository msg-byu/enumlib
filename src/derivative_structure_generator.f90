!!<summary>Generate all derivative structures of a parent structure
!!Gus L. W. Hart BYU July 2007</summary>
MODULE derivative_structure_generator
use utilities
use io_utils
use num_types
use enumeration_types
use labeling_related
use combinatorics
use symmetry
use compare_structures
use rational_mathematics
use numerical_utilities
use vector_matrix_utilities, only: determinant, matrix_inverse
use sorting
use classes, only: polya
use arrow_related, only: arrow_concs

implicit none
private
public get_all_HNFs, remove_duplicate_lattices, get_SNF, get_all_2D_HNFs,&
     &  gen_multilattice_derivatives, find_permutation_of_group,&
     get_dvector_permutations, get_rotation_perms_lists, do_rotperms_form_groups,&
     & mixed_radix_counter

CONTAINS

  !!<summary>This routine builds a concentration table with numerators
  !!and denominators for the given cell size subject to the
  !!concentration constraints specified in cTable. GH Apr 2011 --
  !!extended the routine for multilattices</summary>
  !!<parameter name="vol" regular="true">Current size of cells in the
  !!enumeration loops.</parameter>
  !!<parameter name="nD" regular="true">Number of
  !!d-vectors.</parameter>
  !!<parameter name="concTable" regular="true">Table of desired
  !!concentration ranges</parameter>
  !!<parameter name="cList" regular="true">List of fractions that are
  !!within the ranges.</parameter>  
  !!<parameter name="eps" regular="true">Finite precision
  !!tolerance.</parameter>
  SUBROUTINE get_list_at_this_size(vol,nD,concTable,cList,eps)
    integer, intent(in) :: vol, nD 
    integer, intent(in) :: concTable(:,:) 
    integer, pointer    :: cList(:,:) ! INTENT(OUT): 
    real(dp), intent(in):: eps
    
    integer volTable(size(concTable,1),size(concTable,2))
    integer i, denom, minv, maxv
    
    do i = 1, size(concTable,1) ! Loop over the rows in the
                                ! table. There are k rows (number of
                                ! "colors")
       denom = concTable(i,3) ! Cell size of input concentration ranges
       minv = minval(concTable(i,1:2)); maxv = maxval(concTable(i,1:2))
       ! Use floor on the bottom of the range and ceiling at the top
       ! so that the resulting integers must necessarily straddle the
       ! desired concentration range. Concentrations outside the range
       ! will be eliminated by get_concentration_list
       volTable(i,1:2) = (/ floor(real(minv,dp)/denom*vol*nD), ceiling(real(maxv,dp)/denom*vol*nD) /) 
       volTable(i,3) = vol*nD ! Set the denominator  of the fractions listed
       ! in the volume table to be the current size of unit cells (i.e.,
       ! the volume) being enumerated.
    enddo

    call get_concentration_list(concTable,volTable,cList,eps)

  END SUBROUTINE get_list_at_this_size

  !!<summary>This routine takes the concentration "table" from lat.in
  !!or struct_enum.in and generates a list of concentration vectors
  !!that are consistent with the ranges in the input file.</summary>
  !!<parameter name="concTable" regular="true">concentration ranges
  !!specified in the input (original)</parameter>
  !!<parameter name="volTable" regular="true">concentration ranges
  !!specified for current unit cell sizes.</parameter>
  !!<parameter name="cList">List of concentration vectors that are
  !!within the ranges.</parameter>
  !!<parameter name="eps" regular="true">Finite precision
  !!tolerance.</parameter>
  SUBROUTINE get_concentration_list(concTable,volTable,cList,eps)
    integer, intent(in) :: concTable(:,:) 
    integer, intent(in) :: volTable(:,:) 
    integer, pointer    :: cList(:,:)  
    real(dp), intent(in):: eps

    integer j, k, i, n, cc
    integer, dimension(size(volTable,1)) :: digCnt, digit
    integer, allocatable :: label(:,:), a(:)
    real(dp), dimension(size(volTable,1)) :: minv, maxv, conc
    
    open(98,file="debug_conc_check.out",position="append")
    write(98,'(" Original concentration ranges:")')
    do i = 1, size(concTable,1)
       write(98,'(90(i4,1x))') concTable(i,:)
    enddo
    write(98,'(/," Numbers for given cell size:")')
    do i = 1, size(volTable,1)
       write(98,'(90(i4,1x))') volTable(i,:)
    enddo
    
    ! Imagine an odometer with n wheels, each which has a sequence of
    ! labels (the allowed labels for that atomic site). We'll turn the
    ! wheels on this odometer and keep a list of those reading that
    ! are consistent with (don't violate) the concentration
    ! restrictions.
    n = volTable(1,3) ! Total number of slots in the labeling
    digit = volTable(:,2) - volTable(:,1) + 1 ! Number of "labels" on
    ! each "wheel" of the odometer
    digCnt = 1 ! Start each wheel of the odometer at the first position
    k = size(volTable,1) ! Number of labels, i.e., number of wheels on
    ! the odometer (i.e, the number of atom types)
    allocate(cList(product(digit),k))
    allocate(a(k))  ! The reading on the "odometer"

    ! Define a table that stores the possible labels for each wheel
    allocate(label(size(volTable,1),maxval(digit)))
    label = -1 ! Ends of the rows in this ragged list Paint the labels
    ! on the wheels (leave -1 for wheels that have more positions than
    ! valid labels for that site)
    do j = 1, k 
       label(j,:) = (/(i,i=volTable(j,1),volTable(j,2))/)
       minv(j) = real(minval((/concTable(j,1),concTable(j,2)/)),dp)/concTable(j,3)
       maxv(j) = real(maxval((/concTable(j,1),concTable(j,2)/)),dp)/concTable(j,3)
    enddo
    ! Initialize the reading, "a", to the beginning
    forall(j=1:k);a(j)=label(j,1);endforall

    cc = 0 ! number of valid concentrations
    do
       !write(98,'(10(i1,1x))') a
       if (sum(a)==n) then ! This "reading" is a valid partition, i.e., cell size is correct
          conc = a/real(n,dp)
          if (.not. ( any(conc<minv-eps) .or. any(conc>maxv+eps) )) then 
             cc = cc + 1 ! Keep track of the number of valid concentration vectors
             cList(cc,:) = a ! Store this in the list
          endif
       endif
       j = k  ! Start at the right end of the digit list, j is the digit pointer
       do
          if (digCnt(j) /= digit(j)) exit ! Is this wheel ready to roll over? No then exit
          a(j) = label(j,1) ! Reset this wheel to the first label
          digCnt(j) = 1     ! Reset the counter for this wheel
          j = j - 1         ! Look to the next wheel to the left
          if (j < 1) exit   ! If the whole odometer is ready to roll over the we are done
       enddo
       if (j < 1) exit ! All the numbers have been visited
       digCnt(j) = digCnt(j) + 1 ! Advance this wheel
       a(j) = label(j,digCnt(j)) ! Put the appropriate label in the j-th spot
    enddo
    cList => ralloc(cList,cc,k) ! Reallocate the list to be the proper size
    write(98,'(/," Number of discrete partitions: ",i5)') size(cList,1)
    write(98,'(/," generated list:")')
    do i = 1, size(cList,1)
       write(98,'(10(i4,1x))') cList(i,:)
    enddo
    close(98)
    
  END SUBROUTINE get_concentration_list

  !!<summary>Just a test routine for expanding the code to treat
  !!different label sets on different sites. This is just a counter of
  !!mixed radix form.</summary>
  !!<parameter name="labels" regular="true">Second index is the digit
  !!#, first is the label #.</parameter>
  !!<parameter name="digit" regular="true">List of the upper bound of
  !!each digit.</parameter>
  SUBROUTINE mixed_radix_counter(labels,digit)
    integer :: labels(:,:) 
    integer :: digit(:) 
    integer :: a(size(digit)), counter(size(digit)) ! The odometer; ordinal counter for each digit
    integer :: i,ic,j,n,ilab, quot
    integer :: nUq ! number of unique numbers possible in this mixed radix system
    integer :: multiplier(size(digit)), digIdx(size(digit)), temp(1), b(size(digit)), idx

    n = size(digit) ! Number of sites
    counter = 1

    ! Find the total number of numbers that the mixed radix system can represent
    nUq = product(digit) !;print *," Number of mixed-radix numbers",nUq

    multiplier = 0; multiplier(n)=1
    do i = n-1,1,-1
       multiplier(i) = digit(i+1)*multiplier(i+1)
    enddo
    write(*,'("Multiplier: ",20i6)') multiplier
    
    a = labels(1,:)
    write(*,'("Starting value: ",10i2)') a
    counter(n) = 0; ic = 0
    do; ic = ic + 1
       if(ic > nUq + 1) stop "Fail safe triggered. Wrong number of iterations in mixed-radix counterem"
       j = n ! Reset the digit index; start with j = the last place in the number
       do !;print *, "inner loop"
          if (counter(j) /= digit(j)) exit ! The digit isn't ready to
          ! roll over, so exit and advance it
          counter(j) = 1     ! Roll back j-th digit to the beginning
          j = j - 1          ! Move to the next (on the left) digit
          if (j < 1) exit    ! Leftmost digit updated, we've finished all the numbers
       enddo
       if (j < 1) exit  ! We're done counting, exit
       counter(j) = counter(j) + 1
       forall(i=1:n); a(i) = labels(counter(i),i);endforall
       do i=1,n; temp = minloc(labels(:,i),a(i)==labels(:,i)) - 1; digIdx(i) = temp(1) ;enddo
       idx = sum(digIdx*multiplier)
       idx = sum((counter-1)*multiplier)
       write(*,'(i6,":",2x,10i2)',advance="no") ic,counter
       write(*,'(i10,1x)',advance="no") idx 
       write(*,'(3x,10i2)',advance="no") a
       write(*,'(3x,10i2)',advance="no") digIdx
       do ilab = 1, n ! Loop over each place (digit) in the labeling
          quot = idx/multiplier(ilab) ! divide the index by next highest multiplier
          b(ilab) = labels(quot+1,ilab) ! store the quotient (indexed to labels)
          idx = idx - quot*multiplier(ilab) ! find the remainder
       enddo
       write(*,'(3x,10i2)') b
    enddo
  END SUBROUTINE mixed_radix_counter

  !!<summary>This function checks that every "product" of two
  !!permutations in a list is still in the list. That is that every
  !!list forms a group. For finite groups, this is a sufficient
  !!condition.</summary>
  !!<parameter name="rpl" regular="true">The reduced (unique) set of
  !!rotation permutations lists for a given index</parameter>
  FUNCTION do_rotperms_form_groups(rpl)
    logical do_rotperms_form_groups
    type(RotPermList) :: rpl(:) 

    integer nL, iL, iP, jP, nP, kP, ng
    logical exists
    integer, allocatable :: testperm(:)
    
    nL = size(rpl) ! Number of lists
    ng = size(rpl(1)%perm,2) ! Number of elements in each permutations
    allocate(testperm(ng))
    
    do_rotperms_form_groups = .true.
    lists: do iL = 1, nL ! Check each list in the set
       nP = size(rpl(iL)%perm,1) ! Number of permutations in this list
       write(11,'(3x)',advance="no")
       do iP = 1, nP; write(11,'(i3)',advance="no") iP; enddo; write(11,*) ! Make
       ! the column headings
       perms:   do iP = 1, nP
          write(11,'(i3)',advance="no") iP   ! Write the row heading
          !do jP = 1, nP; write(11,'(3x)',advance='no');enddo ! skip to the column
          do jP = 1, nP ! Now loop over all products for the iP'th permutation
             exists = .false. 
             ! Is the product of perm_i x perm_j in the set?  Permute
             ! the elements of the iP'th permutation according to the
             ! jP'th permutation
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
       if (exists)write(11,*) iL,"The permutations form a group"
       if(.not. exists) write(11,*) iL,"failed"
    enddo lists ! Loop over lists
    !if(exists) then; do_rotperms_form_groups = .true.
    !else; do_rotperms_form_groups = .false.
    !endif
  ENDFUNCTION do_rotperms_form_groups

  !!<summary>This routine takes a list of permutation lists and
  !!identifies those that are idential. The output, RPLindx, is an
  !!index that groups the lists that match. The lists themselves are
  !!assumed to be in "alphabetical" order so that they can be quickly
  !!compared. I don't think it's necessary to explicitly treat lists
  !!corresponding to different SNFs separately---if the permutations
  !!lists happen to be identical (is this possible?) then the
  !!labelings list will also be, so it's not necessary to create a
  !!separate list or index for them.</summary>
  !!<parameter name="RPList" regular="true">The site permutation list.</parameter>
  !!<parameter name="rdRPList" regular="true">The reduced permutation list</parameter>
  !!<parameter name="RPLindx" regular="true"></parameter>
  !!<parameter name="aperms" regular="true">The arrow permutation list.</parameter>
  !!<parameter name="rdaperms" regular="true">The reduced arrow permutation list.</parameter>
  SUBROUTINE organize_rotperm_lists(RPList,rdRPList,RPLindx,aperms,rdaperms)
    type(RotPermList), intent(in) :: RPList(:)
    type(RotPermList), intent(in) ::aperms(:)
    type(RotPermList), pointer, intent(out) :: rdRPList(:), rdaperms(:) ! output
    integer, pointer, intent(out) :: RPLindx(:) ! output

    integer iL, jL, nL, status,  cnt, nP
    !!<local name="taList">Temporary arrow permutation list.</local>
    type(RotPermList), allocatable :: tList(:),taList(:)
    logical unique

    nL = size(RPlist) ! Number of lists (including duplicates)
    allocate(tList(nL),taList(nL),RPLindx(nL),STAT=status)
    if(status/=0) stop "Allocation failed in organize_rotperm_lists: tList, RPLindx"

    cnt = 0
    do iL = 1, nL ! Loop over each loop in the list
       unique = .true.
       do jL = 1, cnt ! Loop over the number of unique lists found so far
          ! if the two lists aren't the same length, then they definitely aren't identical
          if (size(RPList(iL)%perm,1)/=size(tList(jL)%perm,1)) cycle
          if (size(aperms(iL)%perm,1)/=size(taList(jL)%perm,1)) cycle
          if (all(RPList(iL)%perm==tList(jL)%perm) .and. &
               all(aperms(iL)%perm==taList(jL)%perm)) then 
                !identical. Tag it and go to next
                unique = .false.
                exit
          endif
       enddo
       if (unique) then ! store this list in the master list
          cnt = cnt + 1 ! number of unique lists found so far
          nP = size(RPList(iL)%perm,1) ! Number of perms in this list
          allocate(tList(cnt)%perm(nP,size(RPList(iL)%perm,2)))
          allocate(taList(cnt)%perm(nP,size(aperms(1)%perm,2)))
          tList(cnt)%nL = nP ! store number and perms in a temporary
          tList(cnt)%perm = RPlist(iL)%perm ! array
          taList(cnt)%perm = aperms(iL)%perm ! array
       endif
       ! Store the label for this unique list in the output index
       RPLindx(iL) = jL 
    enddo
    ! Now copy the reduced list to the output variable
    allocate(rdRPList(cnt),rdaperms(cnt))
    do iL=1,cnt
       nL = size(tList(iL)%perm,1)
       allocate(rdRPList(iL)%perm(nl,size(tList(iL)%perm,2)));
       allocate(rdaperms(iL)%perm(nl,size(aperms(1)%perm,2)))

       rdRPList(iL)%nL = nL
       rdRPList(iL)%perm = tList(iL)%perm
       if (size(aperms(1)%perm,2) > 1) then
          rdaperms(iL)%perm = taList(iL)%perm
       else
          rdaperms(iL)%perm = 0
       end if
       rdaperms(iL)%nL = nL
       rdaperms(iL)%v => null()
       rdaperms(iL)%RotIndx => null()
       rdRPList(iL)%v => null()
       rdRPList(iL)%RotIndx => null()
    enddo
  ENDSUBROUTINE organize_rotperm_lists
  
  !!<summary>This routine applies the symmetry operations of the
  !!parent lattice to the interior points (i.e., the d-set) to see
  !!which ones are equivalent. Labelings of the lattice points that
  !!are contain permutations only of labels on equivalant sites are
  !!physically equivalent and therefore redundant. We use these
  !!permutations to eliminate those duplicate labelings.</summary>
  !!<parameter name="pLV" regular="true">Lattice vectors of the
  !!primary lattice (parent lattice).</parameter>
  !!<parameter name="pd">d-vectors defining the multilattice (primary
  !!lattice only).</parameter>
  !!<parameter name="dRPList" regular="true">Output. A list
  !!permutations effected by the Ops.</parameter>
  !!<parameter name="Latdim" regular="true">2 or 3 dimensional
  !!case?</parameter>
  !!<parameter name="eps" regular="true">Finite precision
  !!tolerance.</parameter>
  SUBROUTINE get_dvector_permutations(pLV,pd,dRPList,LatDim,eps)
    real(dp) :: pLV(3,3) 
    real(dp), pointer :: pd(:,:) 
    type(RotPermList), intent(out) :: dRPList 
    integer, intent(in) :: LatDim 
    real(dp), intent(in) :: eps 

    integer nD, iD, nOp, iOp, status
    integer, pointer :: aTyp(:), tList(:,:)
    real(dp) :: rd(size(pd,1),size(pd,2)), tRD(size(pd,1),size(pd,2))
    real(dp) :: inv_pLV(3,3) ! Inverse of the pLV matrix
    real(dp), pointer:: rot(:,:,:), shift(:,:), tv(:,:,:)
    logical err
    character(80) name
    nD = size(pd,2)
    allocate(aTyp(nD),STAT=status)
    if(status/=0)stop "Allocation failed in get_dvector_permutations: aTyp"
    
    aTyp = 1
    call get_spaceGroup(pLV,aTyp,pd,rot,shift,.false.,eps)
    
    !call write_lattice_symmetry_ops(rot,shift)
    if(latDim==2) call rm_3D_operations(pLV,rot,shift,eps)
    !call write_lattice_symmetry_ops(rot,shift,"2D")
    nOp = size(rot,3)
    allocate(tList(nOp,nD),tv(3,nD,nOp),STAT=status)
    if(status/=0)stop "Allocation failed in get_dvector_permutations: tList"
    allocate(dRPList%perm(nOp,nD),dRPList%v(3,nD,nOp),STAT=status)
    if(status/=0)stop "Allocation failed in get_dvector_permutations: dRPList"
    
    call matrix_inverse(pLV,inv_pLV,err)
    if (err) stop "Bad parent lattice vectors in input to get_dvector_permutations"
    
    dRPList%nL = nOp  ! Number of operations that fix the parent
    dRPList%RotIndx => null()
    !lattice (but may permute the d-vectors)

    do iOp = 1, nOp ! Try each operation in turn and see how the
                    ! d-vectors are permuted for each
       rd = matmul(rot(:,:,iOp),pd)+spread(shift(:,iOp),2,nD) ! Rotate each d and add the shift
       tRD = rd
       do iD = 1, nD
          call bring_into_cell(rd(:,iD),inv_pLV,pLV,eps)
       enddo
 
       ! The v vector is the vector that must be added (it's a lattice
       ! vector) to move a rotated d-vector back into the parent cell.
       dRPList%v(:,:,iOp) = rd(:,:) - tRD(:,:)
       call map_dvector_permutation(rd,pd,dRPList%perm(iOp,:),eps)
    enddo

    name = "debug_dvec_rots.out"
    call write_rotperms_list(dRPList,name)

    ! I don't think we should reduce this list to a unique one. Some
    ! rotations that don't permute the d's could still permute the
    ! g's. So we have to keep all the d permutations, even if they
    ! look redundant here.

  ENDSUBROUTINE get_dvector_permutations

  !!<summary>For each HNF, we have a list of the operations (rotations
  !!+ shifts, if present) that leave the superlattice fixed. Given
  !!this set of fixing operations, make a list of the permutations on
  !!the d-vectors (interior points of the multilattice) effected by
  !!the rotations. Then sort the HNFs into blocks that have the same
  !!rotation permutations.</summary>
  !!<parameter name="A" regular="true">Lattice vectors of the primary
  !!lattice (parent lattice).</parameter>
  !!<parameter name="HNF" regular="true">List of HNF
  !!matrices.</parameter>
  !!<parameter name="L" regular="true">List of left
  !!transforms.</parameter>
  !!<parameter name="SNF" regular="true">List of SNF
  !!matrices.</parameter>
  !!<parameter name="Op" regular="true">A list of symmetry ops (rots
  !!and shifts) for the parent multilattice</parameter>
  !!<parameter name="RPlist" regular="true">A list of lists of
  !!permutations effected by the Ops.</parameter>
  !!<parameter name="dperms" regular="true"></parameter>
  !!<parameter name="eps" regular="true">Finite precision
  !!tolerance.</parameter>
  !!<parameter name="aperms" regular="true">The output arrow
  !!permutations.</parameter>
  !!<parameter name="use_arrows" regular="true">True if arrows are present
  !!in the enumeration.</parameter>
  !!<parameter name="surf" regular="true">True is this is a surface
  !!calculation, false otherwise.</parameter>
  SUBROUTINE get_rotation_perms_lists(A,HNF,L,SNF,Op,RPlist,dperms,eps,aperms,use_arrows, surf)
    real(dp), intent(in) :: A(3,3) 
    integer, intent(in), dimension(:,:,:) :: HNF, L, SNF 
    type(OpList), intent(in) :: Op(:) 
    type(RotPermList) :: RPlist(:) 
    type(RotPermList), intent(in) :: dperms
    real(dp), intent(in) :: eps
    type(RotPermList), optional, pointer :: aperms(:)
    logical, optional, intent(in) :: use_arrows, surf
    
    !!<local name="taperms">Temporary storage of the arrow permutations.</local>
    !!<local name="arrows">Logigal, true if arrows are to be used.</local>
    type(RotPermList) :: rperms, tperms, taperms
    integer, pointer :: g(:,:) => null()
    integer, allocatable :: gp(:,:), dgp(:,:), ag(:,:), dap(:) ! G prime; the "rotated"
    ! group, (d',g') "rotated" table
    integer, allocatable :: tg(:,:), perm(:), ident(:,:), identT(:,:) ! translated
    ! group, translation permutation of the group members
    integer iH, nH, diag(3), iD, nD, iOp, nOp, n, ig, it, status, im, jm
    real(dp), dimension(3,3) :: Ainv, T, Tinv
    logical err
    real(dp), allocatable :: rgp(:,:), rag(:,:)
    logical, allocatable :: skip(:), skipa(:)
    integer OpIndxInSuperCellList, RowInDxGTable
    integer, allocatable :: arrow_basis(:,:)
    integer, pointer :: temp_perms(:,:) => null()
    logical :: arrows

    if (.not. present(use_arrows)) then
       arrows = .false.
    else
       arrows = use_arrows
    end if

    ! build the arrow basis.
    if (present(surf) .and. (surf .eqv. .True.)) then
       allocate(arrow_basis(3,4))
       arrow_basis(:,1) = (/1,0,0/)
       arrow_basis(:,2) = (/-1,0,0/)
       arrow_basis(:,3) = (/0,1,0/)
       arrow_basis(:,4) = (/0,-1,0/)
    else
       allocate(arrow_basis(3,6))
       arrow_basis(:,1) = (/1,0,0/)
       arrow_basis(:,2) = (/-1,0,0/)
       arrow_basis(:,3) = (/0,1,0/)
       arrow_basis(:,4) = (/0,-1,0/)
       arrow_basis(:,5) = (/0,0,1/)
       arrow_basis(:,6) = (/0,0,-1/)
    end if

    open(19,file="debug_get_rotation_perms_lists.out")

    ! Number of HNFs (superlattices); Index of the superlattices;
    ! Number of d-vectors in d set
    nH = size(HNF,3); n = determinant(HNF(:,:,1)); nD = size(RPList(1)%v,2)

    allocate(gp(3,n), dgp(nD,n), rgp(3,n), ag(3,size(arrow_basis,2)), rag(3,size(arrow_basis,2)), skip(n),dap(size(arrow_basis,2)),skipa(size(arrow_basis,2)),STAT=status)
    if(status/=0) stop "Allocation failed in get_rotation_perm_lists: gp, dgp, rgp, skip" 
    allocate(tg(3,n),perm(n),ident(nD,n),identT(n,nD),STAT=status)
    if(status/=0) stop "Allocation failed in get_rotation_perm_lists: tg, perm, ident"
    if (present(aperms)) then
       allocate(aperms(size(RPlist)),STAT=status)
       ! if(status/=0) then "Allocation failed in get_rotation_perms_lists: aperms."
    end if
    allocate(tperms%perm(n,n*nD),STAT=status)
    ! if(status/=0) then "Allocation failed in get_rotation_perms_lists: tperms."
    identT = reshape((/(ig,ig=1,n*nD)/),(/n,nD/))  ! we could combine
    ! those two lines, but gfortran
    ident  = transpose(identT)                     ! does some strange things then.
    forall(iH = 1:nH); RPlist(iH)%nL=0; end forall ! initialize the number in each list

    ! Make the group member list for the first SNF
    diag = (/SNF(1,1,1),SNF(2,2,1),SNF(3,3,1)/)
    call make_member_list(diag,g)
    call matrix_inverse(A,Ainv,err)
    if(err) stop "Invalid parent lattice vectors in get_rotation_perm_lists"

    do iH = 1,nH ! loop over each superlattice unless the SNF is
   ! different than the previous (they should be sorted into blocks)
   ! don't bother making the group again. Just use the same one.
       if(iH > 1) then; if(.not. all(SNF(:,:,ih)==SNF(:,:,ih-1))) then
          diag = (/SNF(1,1,ih),SNF(2,2,ih),SNF(3,3,ih)/)
          call make_member_list(diag,g)
       endif; endif
       ! Make the transform matrices for taking the g's and rotating them
       Tinv = matmul(L(:,:,iH),Ainv); call matrix_inverse(Tinv, T, err)
       if (err) stop "Bad inverse for transformation matrix: get_rotation_perm_lists"
       
       nOp = size(Op(iH)%rot,3);
       if (associated(rperms%perm)) deallocate(rperms%perm)
       if (associated(taperms%perm)) deallocate(taperms%perm)

       allocate(rperms%perm(nOp,nD*n),taperms%perm(nOp, size(arrow_basis,2)),STAT=status)
       if (status/=0) stop "Allocation failed in get_rotation_perm_lists: rperms%perm"
       do iOp = 1, nOp ! For each rotation, find the permutation
          dgp = 0 ! Initialize the (d,g) table
          dap = 0
          do iD = 1, nD ! Loop over each row in the (d,g) table
             ! LA^-1(v_i+(RAL^-1)G)
             rgp = matmul(Tinv,(-spread(RPList(iH)%v(:,iD,iOp),2,n)+matmul(matmul(Op(iH)%rot(:,:,iOp),T),g)))
             rag = transpose(matmul(transpose(arrow_basis),op(iH)%rot(:,:,iOp)))
             if (.not. equal(rgp,nint(rgp),eps)) stop "Transform left big fractional parts"
             gp = nint(rgp) ! Move the rotated group into an integer array
             ag = nint(rag)
             gp = modulo(gp,spread(diag,2,n)) ! Mod by each entry of
             ! the SNF to bring into group Now that the rotated group
             ! is known, find the mapping of the elements between the
             ! original group and the permuted group. This is the
             ! permutation.
             skip = .false. ! This is just for efficiency
             do im = 1, n
                do jm = 1, n
                   if (skip(jm)) cycle ! Skip elements whose mapping is already known
                   if (all(gp(:,jm)==g(:,im))) then ! these elements
                   ! map to each other The list of operations that fix
                   ! the superlattice are a subset of those that fix
                   ! the parent lattice. RotIndx stores the indicies
                   ! of the parent lattice operations in a list with
                   ! as many entries as supercell fixing operations.

                   ! dperms%perm stores a list of d-vector
                   ! permutations, one permutation (an nD list) for
                   ! each operation in the parent lattice symmetries
                      OpIndxInSuperCellList = RPList(iH)%RotIndx(iOp) 
                      RowInDxGTable = dperms%perm(OpIndxInSuperCellList,iD)
                      dgp(RowInDxGTable,im) = jm+(iD-1)*n
                      skip(jm) = .true.
                      exit
                   endif
                enddo ! jm
             enddo ! im

             ! do the some operations for the arrows
             skipa = .false.
             do im = 1, size(arrow_basis,2)
                do jm = 1, size(arrow_basis,2)
                   if (skipa(jm)) cycle
                   if (all(ag(:,jm)==arrow_basis(:,im))) then
                      dap(im) = jm
                      skipa(jm) = .true.
                      exit
                   end if
                end do !jm
             end do !im
          enddo ! loop over d-vectors (each row in the table)
          if (any(dgp==0) .or. (any(dap==0) .and. arrows)) stop "(d,g)-->(d',g') mapping failed in get_rotation_perm_lists"

          ! Now we have the (d',g') table for this rotation. Now
          ! record the permutation
          rperms%perm(iOp,:) = reshape(transpose(dgp),(/nD*n/)) ! store
          taperms%perm(iOp,:) = dap
          ! permutation in the "long form"
          rperms%nL = nOp
          !write(*,'(i2,1x,20i2)')iOp,rperms%perm(iOp,:)
          !write(*,'("nL",1x,20i2)')rperms%nL
       enddo ! loop over rotations 

       ! nomenclature:
       ! N+t = rotation (N) + fractional translation (t)  (me bethinks....)
       !   r = lattice translation

       ! Now that we have the permutations that are effected by N+t
       ! type of rotations (for each N in the point group), we need to
       ! compose them with the permutations effected by lattice
       ! translations, (parent) lattice translations contained inside the
       ! supercell (N+t+r, in Rod's nomenclature). Only when we have these
       ! two lists, and compose them, do we have a guarantee that the
       ! resulting list is actually a group. For efficiency, we reduce the
       ! N+t list (remove duplicates). Rod claims that the compositions
       ! (with the translations) will not have duplicates.
       if ((size(rperms%perm,1) > 1) .and. (arrows .eqv. .false.)) then
          call sort_permutations_list(rperms%perm)

       else if (arrows .eqv. .True.) then
          allocate(temp_perms(size(rperms%perm,1),size(rperms%perm,2)+size(taperms%perm,2)))
          do im=1, size(rperms%perm,1)
             temp_perms(im,1:size(rperms%perm,2)) = rperms%perm(im,:)
             temp_perms(im,size(rperms%perm,2)+1:) = taperms%perm(im,:)
          end do
          call sort_permutations_list(temp_perms)
          deallocate(rperms%perm,taperms%perm)
          allocate(rperms%perm(size(temp_perms,1),size(temp_perms,2)-size(arrow_basis,2)),STAT=status)
          if(status/=0) stop "Allocation failed in get_rotation_perms_lists: rperms"
          allocate(taperms%perm(size(temp_perms,1),size(arrow_basis,2)),STAT=status)
          if(status/=0) stop "Allocation failed in get_rotation_perms_lists: taperms"
          do im =1,size(temp_perms,1)
             rperms%perm(im,:) = temp_perms(im,1:size(rperms%perm,2))
             taperms%perm(im,:) = temp_perms(im,size(rperms%perm,2)+1:)
          end do
       end if
       ! The rotations permutations list is now in "alphabetical"
       ! order and contains no duplicates

       ! To get the permutations effected by the r's, we don't need
       ! the r's. We can merely take the member list, the group, and
       ! add each element of the group to the group itself and see
       ! what permutation happens. (This part is somewhat redundant
       ! with make_translation_group in labeling_related module.)
       do ig = 1,n ! The number of r's inside the superlattice (the
                   ! translation perms) is the same as the index n
          tg = g(:,:)+spread(g(:,ig),2,n) ! Add the element to the group
          tg = modulo(tg,spread(diag,2,n)) ! mod by the SNF entries to
          ! bring it back to the "primitive" representation
          call find_permutation_of_group(g,tg,perm)
          tperms%perm(ig,:) = reshape(transpose(ident(:,perm)),(/n*nD/))
       enddo
       
       RPlist(iH)%nL = size(rperms%perm,1)*n
       allocate(RPlist(iH)%perm(RPlist(iH)%nL,n*nD)) ! nL rows and n*nD columns in the list
       if ((arrows) .and. present(aperms)) then
          allocate(aperms(iH)%perm(RPlist(iH)%nL,size(arrow_basis,2)))
       else if (present(aperms)) then
          allocate(aperms(iH)%perm(RPlist(iH)%nL,1))
       end if
       do it = 1,n ! Loop over translation perms (r type)
          do iOp = 1,size(rperms%perm,1) ! Loop over unique rotation
             ! perms (N+t type) Form the permutation effected by
             ! composing the iOp-th one with the it-th one
             RPlist(iH)%perm((iOp-1)*n+it,:) = tperms%perm(it,(rperms%perm(iOp,:)))
             if ((arrows .eqv. .true.) .and. present(aperms)) then
                aperms(iH)%perm((iOp-1)*n+it,:) = taperms%perm(iOp,:)
             else if (present(aperms)) then
                aperms(iH)%perm((iOp-1)*n+it,:) = 0
             end if
             ! ^--- Having gotten both the rotations and the
         ! translation in the sections above (sort_permutations_list
         ! etc...), the "operators" of the rotation (one of them is
         ! R_i) and the translations (one of them is T_j) are now
         ! "scrambled", i.e., (T_j) o (R_i).  Since the *first*
         ! Rotation R_1 = Id, the *first* entries in the
         ! RPlist(iH)%perm are equivalent to pure translations only.
         ! The following entries are a combination of both.
          enddo
       enddo
    enddo ! loop over iH (superlattices)
    close(19)
  ENDSUBROUTINE get_rotation_perms_lists

  !!<summary></summary>
  !!<parameter name="g" regular="true">Unpermuted groups.</parameter>
  !!<parameter name="gp" regular="true">Permuted groups.</parameter>
  !!<parameter name="perm" regular="true">Permutation of gp.</parameter>
  SUBROUTINE find_permutation_of_group(g,gp,perm)
    integer, intent(in), dimension(:,:) :: g, gp 
    integer, intent(out) :: perm(:) 

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

  !!<summary></summary>
  !!<parameter name="rd" regular="true">Rotated d's</parameter>
  !!<parameter name="d" regular="true">Original d.</parameter>
  !!<parameter name="RP" regular="true">Permutation of the quotient
  !!group effected by the rotation.</parameter>
  !!<parameter name="eps" regular="true">Finite precision tolerance.</parameter>
  SUBROUTINE map_dvector_permutation(rd,d,RP,eps)
    real(dp), dimension(:,:), intent(in) :: rd, d ! (both input)
    integer, intent(out) :: RP(:) !(output)
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
    if(any(RP==0)) then; print *, "d-vector didn't permute in map_dvector_permutation";
       print *,"This usually means that the d-set from the input structure and the d-set"
       print *,"from the struct_enum.out have a different origin or don't live in the same"
       print *,"unit cell. This probably isn't your fault---the code should overcome this."
       write(*,'(200(i2,1x))') RP
    stop;endif
  ENDSUBROUTINE map_dvector_permutation

  !!<summary>Finds all the possible diagonals of the HNF matrices of a
  !!given size</summary>
  !!<parameter name="detS" regular="true">Cell size, i.e., determinant
  !!of S matrix.</parameter>
  !!<parameter name="diagonals">All possible diagonals.</parameter>
  SUBROUTINE get_HNF_diagonals(detS, diagonals)
    integer, intent(in) :: detS 
    integer, pointer :: diagonals(:,:) 

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
  
  !!<summary>This subroutine generates all the unique HNF matrices of
  !!a given determinant. (See Santoro and Mighell 1972
  !!Acta. Cryst.)</summary>
  !!<parameter name="volume" regular="true"></parameter>
  !!<parameter name="hnf"></parameter>
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
                hnf(:,:,ihnf) = reshape((/ d(1,i),      j,     k,        &   
                     0, d(2,i),     l,        &   
                     0,      0, d(3,i)  /), (/3,3/))
       enddo;enddo;enddo  ! End loops over values for off-diagonal elements
    enddo ! End loop over all unique triplets of target determinant (volume)

    if (ihnf /= Nhnf) stop "HNF: not all the matrices were generated...(bug!)"
  END SUBROUTINE get_all_HNFs

  !!<summary>Find all the SNFs of a group of HNFs. A and B are the
  !!transformation matrices. F is a list of the *unique*
  !!SNFs. SNF_label indicates which of the unique SNFs corresponds to
  !!each HNF. The transformations, SNFs, and labels are ordered on
  !!output into blocks of identical SNFs.</summary>
  !!<parameter name="HNF">List of HNFs (input).</parameter>
  !!<parameter name="A">Transformation matrix (output)</parameter>
  !!<parameter name="SNF">List of SNFs</parameter>
  !!<parameter name="B">Transformation matrix (ouput)</parameter>
  !!<parameter name="RPList" regular="true">List of rotation
  !!permutations for each HNF.</parameter>
  !!<parameter name="F"></parameter>
  !!<parameter name="SNF_label"></parameter>
  !!<parameter name="fixing_op" regular="true">List of operations that
  !!fix each HNF.</parameter>
  SUBROUTINE get_SNF(HNF,A,SNF,B,RPList,F,SNF_label,fixing_op)
    integer, pointer :: HNF(:,:,:), SNF_label(:) ! HNF is an input list
    integer, pointer, dimension(:,:,:) :: A, SNF, B, F ! All these are output lists
    type(opList), intent(inout) :: fixing_op(:) 
    type(RotPermList), intent(inout) :: RPList(:) 

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
    
    ! We have all the SNFs now, as well as the transformations. But
    ! they are not in any kind of order. Rearrange the
    ! transformations, labels, and HNFs so that they are in blocks of
    ! matching SNFs
    ic = 1; jc = 0
    do ifound = 1, nfound
       jc = count(SNF_label==ifound) ! How many HNFs have this SNF?
       indx(ic:ic+jc-1) = pack((/(i,i=1,nHNF)/),SNF_label==ifound) ! Get
       ! the index of each matching case
       ic = ic + jc ! Compute the offset where the index is stored
    enddo

    ! Use the computed order index, indx, to reorder the A, B, and SNF matrices
    HNF = HNF(1:3,1:3,indx)
    SNF = SNF(1:3,1:3,indx)
    A = A(1:3,1:3,indx)
    B = B(1:3,1:3,indx)
    SNF_label = SNF_label(indx) ! Reorder the labels for which SNF each HNF is
    fixing_op = fixing_op(indx)
    RPList = RPList(indx) ! Reorders the permutations lists
    
    ! Fail safe trigger and storage of unique SNFs
    if (ic/=nHNF+1) stop "SNF sort in get_SNF didn't work"
    !if (associated(F)) deallocate(F). This is unnecessary because F
    !should only be allocated here and it only happens once in each
    !call to gen_multilattice_derivatives
    allocate(F(3,3,nfound),STAT=status)
    if(status/=0) stop "Failed to allocate memory in get_SNF for array F"
    F = tF(:,:,1:nfound)
    
  ENDSUBROUTINE get_SNF

  !!<summary>Takes a bunch of HNF matrices and a parent lattice and
  !!removes those that are rotationally equivalent (under the
  !!rotations of the parent lattice). Also returns all the unique
  !!derivative _lattices_ for this parent. Also returns a list of
  !!rotations that fix each superlattice.</summary>
  !!<parameter name="hnf">HNF matrices.</parameter>
  !!<parameter name="LatDim" regular="true">Is the parent lattice 2D
  !!or 3D?</parameter>
  !!<parameter name="parent_lattice" regular="true">Parent
  !!lattice.</parameter>
  !!<parameter name="d"></parameter>
  !!<parameter name="dperms" regular="true"></parameter>
  !!<parameter name="uq_hnf">List of symmetrically distinct
  !!HNFs.</parameter>
  !!<parameter name="fixing_op" regular="true">List of operations that
  !!fix each HNF.</parameter>
  !!<parameter name="RPList"></parameter>
  !!<parameter name="latts"></parameter>
  !!<parameter name="degeneracy_list" regular="true"></parameter>
  !!<parameter name="eps" regular="true">Finite precision tolerance.</parameter>
  SUBROUTINE remove_duplicate_lattices(hnf,LatDim,parent_lattice,d,dperms,uq_hnf,fixing_op,RPList,latts,degeneracy_list,eps)
    integer, pointer :: hnf(:,:,:) ! (input)
    integer :: LatDim 
    real(dp), intent(in) :: parent_lattice(3,3) ! parent lattice (input)
    integer, pointer     :: degeneracy_list(:)
    integer, pointer :: uq_hnf(:,:,:) ! (output)
    type(opList), pointer :: fixing_op(:) 
    real(dp),intent(in):: eps ! finite precision (input)
    real(dp), pointer :: d(:,:) 
    type(RotPermList), intent(in) :: dperms
    type(RotPermList), pointer :: RPList(:)
    real(dp), pointer :: latts(:,:,:)
    !integer, intent(in)  :: label(:,:), digit(:)
    
    real(dp), pointer:: sgrots(:,:,:), sgshift(:,:)
    real(dp), dimension(3,3) :: test_latticei, test_latticej
    integer i, Nhnf, iuq, irot, j, nRot, Nq, status, nD
    integer, allocatable :: temp_hnf(:,:,:)
    logical duplicate
    integer, pointer :: aTyp(:)
    
    nD = size(d,2)
    allocate(aTyp(nD),STAT=status)
    if(status/=0)stop "Allocation failed in remove_duplicate_lattices: aTyp"
    
    
    aTyp = 1
    ! Let the code apply all the symmetries and then eliminate invalid
    !  labelings (the counter won't generate any but they may appear
    !  after the symmetry has been applied to a legal one). In other
    !  words, don't try and restrict the set of symmetries (at this
    !  point) but the set of labels that can be applied to any
    !  particular site.  call
    !  get_spaceGroup_atomTypes(label,digit,aTyp)
    call get_spaceGroup(parent_lattice,aTyp,d,sgrots,sgshift,.false.,eps)
    nRot = size(sgrots,3)
    Nhnf = size(hnf,3)
    allocate(temp_hnf(3,3,Nhnf),STAT=status)
    
    
    if(status/=0) stop "Failed to allocate memory in remove_duplicate_lattices: temp_hnf"
    temp_hnf = hnf
    call write_lattice_symmetry_ops(sgrots,sgshift)
    
    ! for the 2D case, eliminate the "3D" operations.
   
    if (LatDim==2) then
       call rm_3d_operations(parent_lattice,sgrots,sgshift,eps)
       nRot = size(sgrots,3)
       call write_lattice_symmetry_ops(sgrots,sgshift,"2D")
    endif
    ! For each HNF in the list, see if it is a derivative lattice of a
    ! preceding HNF in the list. If so, don't include it in the list
    ! of unique ones.
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
          temp_hnf(:,:,iuq) = hnf(:,:,i)
       endif
      
    enddo

    allocate(uq_hnf(3,3,iuq),latts(3,3,iuq),RPList(iuq),STAT=status)
    if(status/=0) stop "Failed to allocate uq_hnf/latts in remove_duplicate_lattices: uq_hnf"
    uq_hnf = temp_hnf(:,:,:iuq)
    Nq = iuq
    forall (i=1:Nq); latts(:,:,i) = matmul(parent_lattice,uq_hnf(:,:,i));end forall
    ! Now loop through the list of unique superlattices and record
    ! which of the parent lattice symmetry operations leave the
    ! superlattice unchanged. These operations, while they leave the
    ! superlattice unchanged can permute the labels inside. Thus, this
    ! is another source of duplicates.
    allocate(fixing_op(Nq),STAT=status)
    if(status/=0) stop "Allocation of fixing_op failed in remove_duplicate_lattices"
    
    allocate(degeneracy_list(Nq) )
    degeneracy_list = 0
    do iuq = 1, Nq  ! Loop over each unique HNF Determine which
       ! operations in the sym ops of the parent lattice leave the
       ! superlattice fixed. These operations may permute the labeling
       ! so we need them when finding duplicate labelings
       call get_sLV_fixing_operations(uq_hnf(:,:,iuq),parent_lattice,nD,sgrots,sgshift,&
        dperms,fixing_op(iuq),RPList(iuq),degeneracy_list(iuq),eps)
    enddo
    degeneracy_list = degeneracy_list + 1
  END SUBROUTINE remove_duplicate_lattices

  !!<summary></summary>
  !!<parameter name="HNF" regular="true"></parameter>
  !!<parameter name="pLV" regular="true"></parameter>
  !!<parameter name="nD" regular="true"></parameter>
  !!<parameter name="rot" regular="true"></parameter>
  !!<parameter name="shift" regular="true"></parameter>
  !!<parameter name="dPerm" regular="true"></parameter>
  !!<parameter name="fixOp" regular="true">Stores the operations that
  !!leave sLV fixed.</parameter>
  !!<parameter name="rotPerm" regular="true"></parameter>
  !!<parameter name="degeneracy" regular="true"></parameter>
  !!<parameter name="eps" regular="true">Finite precision tolerance.</parameter>
  SUBROUTINE get_sLV_fixing_operations(HNF,pLV,nD,rot,shift,dPerm,fixOp,rotPerm,degeneracy,eps)
    integer, intent(in)       :: HNF(:,:), nD
    real(dp), intent(in)      :: pLV(:,:), rot(:,:,:), shift(:,:)
    type(opList)              :: fixOp 
    real(dp), intent(in)      :: eps                 
    type(RotPermList)         :: rotPerm
    type(RotPermList), intent(in):: dPerm
    
    integer :: degeneracy
    integer i, ic, iRot, nRot, status, iDegen, cDegen
    type(opList)               :: tmpOp
    real(dp), dimension(3,3)   :: thisRot, origLat, rotLat
    real(dp), allocatable      :: tv(:,:,:), degen_lattices(:,:,:)
    integer, allocatable       :: tIndex(:)
    logical :: inList

    nRot = size(rot,3)
    allocate(tv(3,nD,nRot),tIndex(nRot),STAT=status); if(status/=0) stop "tv didn't allocate"
    allocate(tmpOp%rot(3,3,nRot),tmpOp%shift(3,nRot))
    
    allocate(degen_lattices(3,3,nRot) )
    degen_lattices = 0
    cDegen = 0
    ic = 0 ! Counter for the fixing operations
    tv = 0; tIndex = 0 ! temp variables
    do iRot = 1,nRot  ! Loop over each rotation
       thisRot = rot(:,:,iRot) ! Store the rotation
       origLat = matmul(pLV,HNF)  ! Compute the superlattice
       rotLat = matmul(thisRot,origLat)          ! Compute the rotated superlattice
       if (is_equiv_lattice(rotLat,origLat,eps)) then ! this operation
       ! fixes the lattice and should be recorded
          ic = ic + 1
          tmpOp%rot(:,:,ic) = thisRot
          tmpOp%shift(:,ic) = shift(:,iRot)
          tv(:,:,ic) = dPerm%v(:,:,iRot)
          tIndex(ic) = iRot
          ! Added by LN from here
       else
          inList = .false.
          do iDegen = 1, cDegen
             if (is_equiv_lattice(degen_lattices(:,:,iDegen),rotLat,eps) ) then
                inList = .true.
                exit
             end if
          end do
          
          if (.not. inList) then
             cDegen = cDegen + 1
             degen_lattices(:,:,cDegen) = rotLat
          end if
       endif
       ! End of LN additions for degeneracy
    enddo ! Now we know which rotations fix the lattice and how many
          ! there are so store them
    degeneracy = cDegen
    do i=1,ic; allocate(fixOp%rot(3,3,ic),fixOp%shift(3,ic),STAT=status)
       if(status/=0) stop "Allocation of fixing_op(iuq) failed, module deriv..."; enddo
    ! Allocate the storage for them
    fixOp%rot =   tmpOp%rot(:,:,1:ic) ! Stuff the rotations into the permanent array
    fixOp%shift = tmpOp%shift(:,1:ic) ! Stuff the shifts into the permanent array
    allocate(rotPerm%v(3,nD,ic),rotPerm%RotIndx(ic),STAT=status)
    if (status/=0) stop "Allocation of RPList failed in remove_duplicate_lattices" 
    rotPerm%v = tv(:,:,1:ic)
    rotPerm%RotIndx = tIndex(1:ic)
    

  END SUBROUTINE get_sLV_fixing_operations

  !!<summary>This subroutine generates all the unique HNF matrices of
  !!a given determinant (See Santoro and Mighell 1972 Acta. Cryst.)
  !!but for a quasi-2D case</summary>
  !!<parameter name="volume" regular="true"></parameter>
  !!<parameter name="hnf" ></parameter>
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
               0, d(2,i),     j,        &   
               0,      0, d(3,i)  /), (/3,3/))
       enddo  ! End loops over values for off-diagonal elements
    enddo ! End loop over all unique triplets of target determinant (volume)
    if (ihnf /= Nhnf) stop "HNF: not all the matrices were generated...(bug!)"
  END SUBROUTINE get_all_2D_HNFs

  !!<summary>Finds all the possible diagonals of the HNF matrices of a
  !!given size.</summary>
  !!<parameter name="detS" regular="true">Cell size, i.e., determinant
  !!of S matrix.</parameter>
  !!<parameter name="diagonals">All possible diagonals.</parameter>
  SUBROUTINE get_HNF_2D_diagonals(detS, diagonals)
    integer, intent(in) :: detS 
    integer, pointer :: diagonals(:,:) 
    
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

  ![TODO] Add line enumerations for proteins etc.
  !!<summary>This procedure generates derivative superstructures from
  !!multilattices.</summary>
  !!<comments> This is more general and should therefore work on
  !! "mono"-lattices, as the original routine did. The algorithm
  !! implemented here has been slightly reordered from the original
  !! routine discussed in the first paper (Hart and Forcade 2008).
  !! </comments>
  !!<revision date="Oct 2010" author="GH">This routine has long since
  !!replaced the original.</revision>
  !!<revision date="Oct 2010" author="GH">Added the ability to
  !!enumerate for a fixed concentration (EnumIII).</revision>
  !!<parameter name="title" regular="true">A string to identify the
  !!run, echoed in the output file.</parameter>
  !!<parameter name="parLV" regular="true">The parent lattice
  !!vectors.</parameter>
  !!<parameter name="nDFull" regular="true">The number of d-vectors in
  !!the full d-set @CREF[param.dFull].</parameter>
  !!<parameter name="dFull">A list of all the d-vectors including
  !!those with no configurational degrees of freedom.</parameter>
  !!<parameter name="k" regular="true">The number of different kind of
  !!atoms that can be placed on a site (a.k.a. colors).</parameter>
  !!<parameter name="nMin" regular="true">The lower bound on the range
  !!of cell sizes to enumerate over.</parameter>
  !!<parameter name="nMax" regular="true">The upper bound on the range
  !!of cell sizes to enumerate over.</parameter>
  !!<parameter name="pLatTyp" choices="S,B" regular="true"> The parent
  !! lattice type. Possible choices are S = "Surface", B = "Bulk".
  !! </parameter>
  !!<parameter name="eps" regular="true">Tolerance for the finite
  !!difference calculations.</parameter>
  !!<parameter name="full" regular="true"> The full setting returns
  !! all possible combinations of colors at sites even if they are
  !! crystallographically equivalent. Partial (i.e. .false.) removes
  !! combinations of colors where a crystallographically equivalent
  !! site already exists (even though the energy may be
  !! different..</parameter>
  !!<parameter name="labelFull">A list of the labels (index 1) allowed
  !!on each site (index 2).</parameter>
  !!<parameter name="digitFull">A list of the *number* of labels on
  !!each site </parameter>
  !!<parameter name="equivalencies">A listing of sites that need to
  !!have the same color.</parameter>
  !!<parameter name="conc_check" regular="true">Specifies whether the
  !! enumeration should be done in a restricted concentration
  !! range. The restricted range is specified in
  !! @CREF[param.cRange]</parameter>
  !!<parameter name="cRange"> <summary>A restricted concentration
  !!  range for the enumeration. If any concentrations are restricted,
  !!  ranges must be specified for ALL colors (atom types) in the
  !!  enumeration.</summary> <dimension>Rows specify range
  !!  fractions. Row 1: Numerator (START); Row 2: Numerator (END); Row
  !!  3: Denominator.</dimension> <dimension>Columns specify colors
  !!  for range limitations.</dimension> </parameter>
  !!<parameter name="polya" regular="true">Logical flag for if the
  !!code is running the polya mode.</parameter>
 SUBROUTINE gen_multilattice_derivatives(title, parLV, nDFull, dFull, k, nMin, nMax, pLatTyp, eps, full,&
       & labelFull,digitFull,equivalencies, conc_check,cRange, polya)
    integer, intent(in) :: k, nMin, nMax, nDFull
    integer             :: nD
    !Had to change character to 80 from 10 to match the definition in io_utils.read_input
    character(80), intent(in) :: title
    real(dp), intent(in) :: parLV(3,3), eps
    real(dp), pointer :: dFull(:,:), d(:,:)
    character(1), intent(in) :: pLatTyp
    logical, intent(in) :: full 
    integer, intent(inout) :: labelFull(:,:) !
    integer, intent(inout) :: digitFull(:)   !
    integer, allocatable    :: label(:,:), digit(:)
    integer, intent(in) :: equivalencies(:)
    logical, intent(in) :: conc_check
    integer, optional,intent(in)   :: cRange(:,:)
    logical, optional, intent(in) :: polya
    
    integer :: iD, i, ivol, LatDim, Scnt, Tcnt, iBlock, HNFcnt, status, iC, j
    integer, pointer, dimension(:,:,:) :: HNF => null(),SNF => null(), L => null(), R => null()
    integer, pointer :: SNF_labels(:) =>null(), uqSNF(:,:,:) => null()
    integer, pointer, dimension(:,:,:) :: rdHNF =>null()
    real(dp) tstart, tend!, removetime, organizetime, writetime,hnftime,snftime, permtime
    real(dp), allocatable :: tempD(:,:) ! Temporary list of the d-vectors
    real(dp)              :: inv_parLV(3,3) ! Parent lattice inverse matrix
    type(opList), pointer :: fixOp(:)  ! Symmetry operations that leave a multilattice unchanged
    type(RotPermList), pointer :: RPList(:) ! Master list of the rotation permutation lists
    type(RotPermList), pointer :: rdRPList(:) ! Master list of the *unique* rotation
    ! permutation lists
    type(RotPermList) :: ParRPList ! Just the list for the parent lattice
    integer, pointer ::  RPLindx(:) => null() ! Index showing which list of rotation
    !permutations corresponds to which HNF
    real(dp), pointer :: uqlatts(:,:,:) => null()
    character, pointer :: lm(:) ! labeling markers (use to generate the labelings
    !when writing results)
    character(80) filename ! String to pass filenames into output writing routines
    character(200) formatstring
    logical fixed_cells ! This is set to true if we are giving a list of cells from a file
    !instead of generating them all
    integer, pointer :: iRange(:,:) ! Rows: List of the number of atoms of each type
    !("color vector"). Cols: sweep over concentration range.
    logical err
    integer, pointer :: hnf_degen(:), lab_degen(:)
    integer, allocatable :: site_res(:,:)
    !!<local name="arrows">Logical that indicates if arrows are to be
    !!included in the enumeration.</local>
    !!<local name="aperms">The full list of the arrow permutations.</local>
    !!<local name="rdaperms">The reduced list of arrow permutations.</local>
    !!<local name="max_binomial">The largest allowed binomial.</local>
    logical :: arrows
    type(RotPermList), pointer :: aperms(:), rdaperms(:)
    real(dp) :: max_binomial

    max_binomial = 1E10

    ! Divide the dset into members that are enumerated and those that are not
    nD = count( (/(i,i=1,nDFull)/)==equivalencies)
    allocate(d(3,nD), label(size(labelFull,1),nD), digit(nD),tempD(3,nD))
    
    nD = 0
    do iD=1,nDFull
       if (iD==equivalencies(iD)) then ! this dset member is a unique point
          nD=nD+1
          d(:,nD)     = dFull(:,iD)
          label(:,nD) = labelFull(:,iD)
          digit(nD)   = digitFull(iD)
       else ! this dset member is equivalent (concerning the
            ! enumeration!) to a different point.  => Force the label
            ! and digit arrays to be equal for equivalent sites
          labelFull(:,iD) = labelFull(:,equivalencies(iD))
          digitFull(iD)   = digitFull(equivalencies(iD))
       endif
    enddo
    
    ! Are we going to use a set of fixed cells? Or loop over all possible cells?
    fixed_cells = .false.
    open(43,file="fixed_cells.in",status="old",iostat=status)
    if(status==0) then
       write(*,'(A)') "---------------------------------------------------------------------------------------------"
       write(*,'("Generating permutations for fixed cells. index n=",i2," to ",i2)') nMin, nMax
       write(*,'("Be aware: non-primitive structures (super-periodic configurations) are included",/,&
            &"in the final list in this mode.")')
       fixed_cells=.true.
    endif
    close(43)
    
    write(*,'(A)') "---------------------------------------------------------------------------------------------"
    write(*,'("Calculating derivative structures for index n=",i2," to ",i2)') nMin, nMax
    if (conc_check) then
       write(*,'(A)') "Including only structures of which the concentration &
            &of each atom is in the range:"
       do i = 1, k
          write(formatstring,'(A,i1,A,i1,A,i1,A,i1,A)') '("Type:",i2,": ",&
               &i',int(log10(real(cRange(i,1))+1))+1,&
               &',"/",i',int(log10(real(cRange(i,3))+1))+1,',"--",i'&
               &,int(log10(real(cRange(i,2))+1))+1,',"/",i'&
               &,int(log10(real(cRange(i,3))+1))+1,')'
          write(*,formatstring) i,cRange(i,1),cRange(i,3),cRange(i,2),cRange(i,3)
       enddo
    endif
    
    write(*,'("Volume",7x,"CPU",8x,"#HNFs",2x,"#SNFs",&
         &4x,"#reduced",4x,"% dups",6x,"volTot",6x,"RunTot")')
    
    ! Check to see if there are arrows present in the system.
    inquire(FILE="arrows.in",EXIST=arrows)

    ! Set up the output file and write the lattice information
    open(14,file="struct_enum.out")
    !Write the fortpy version information for the file.
    if (arrows) then
       write(14, *) '# <fortpy version="4" revision="247"></fortpy>'
    else
       write(14, *) '# <fortpy version="3" revision="247"></fortpy>'
    end if
    write(14,'(a10)') title
    if (pLatTyp=='S'.or.pLatTyp=="s") then; write(14,'(a4)') "surf"
    elseif (pLatTyp=='B'.or.pLatTyp=="b") then; write(14,'(a4)') "bulk"
    else; print *,"pLatTyp:",pLatTyp
    stop '"pLatTyp" not defined in gen_multilattice_derivs in enumlib'; endif
       
    do i = 1,3
       write(14,'(3(g15.8,1x),3x,"# a",i1," parent lattice vector")') parLV(:,i),i
    enddo
    write(14, '(i5," # Number of points in the multilattice")') nDFull
    do iD = 1,nDFull
       ! Print out the dset points and their possible labels 
       ! (1) setup the format for output      
       formatstring='(3(g15.8,1x),3x,"# d",i2.2," d-vector, labels: "'
       do i=1,digitFull(iD); if (i>1) formatstring=trim(formatstring)//',"/"'; formatstring=trim(formatstring)//',i1'; enddo
       formatstring=trim(formatstring)//")"
       ! (2) print the data
       write(14,formatstring) dFull(:,iD),iD, labelFull(1:digitFull(iD),iD)
    enddo
    write(14,'(i2,"-nary case")') k
    write(14,'(2i4," # Starting and ending cell sizes for search")') nMin, nMax
    write(14,'(g15.8," # Epsilon (finite precision parameter)")') eps
    write(14,'(A)') "Concentration check:"
    write(14,'(L5)') conc_check
    if (conc_check) then
       write(14,'(A)') "Including only structures of which the concentration &
            &  of each atom is in the range:"
       do i = 1, k
          write(14,'("Type ",i1,": ",i4,"/",i4," -- ",i4,"/",i4)') i,cRange(i,1),cRange(i,3),cRange(i,2),cRange(i,3)
       enddo
    endif
    if (full) then; write(14,'("full list of labelings (including incomplete labelings) is used")')
    else; write(14,'("partial list of labelings (complete labelings only) is used")'); endif
    write(14,'("Equivalency list:" ,40(I2,1x))') equivalencies(:)
             
    !write(14,'("Symmetry of the primary lattice is of order ",i2)')

    if (arrows) then
       write(14,'("start",3x,"#tot",6x,"HNF",5x,"Hdegn",3x,"labdegn",3x,"Totdegn",3x,"#size",1x,"idx",4x,"pg",4x,"SNF",13x,"HNF",17x,"Left transform",26x,"labeling",5x,"arrows")')
    else
       write(14,'("start",3x,"#tot",6x,"HNF",5x,"Hdegn",3x,"labdegn",3x,"Totdegn",3x,"#size",1x,"idx",4x,"pg",4x,"SNF",13x,"HNF",17x,"Left transform",26x,"labeling")')
    end if

    ! Check for 2D or 3D request
    if (pLatTyp=='s' .or. pLatTyp=='S') then; LatDim = 2
       if (.not. equal((/parLV(1,2),parLV(1,3)/),(/0._dp,0._dp/),eps)) &
            stop 'For "surf" setting, first component of second and third vectors &
            & must be zero'
    else if(pLatTyp=='b' .or. pLatTyp=='B') then; LatDim = 3
    else; stop 'Specify "surf" or "bulk" in call to "generate_derivative_structures"';endif
                
    ! The permutations of the interior points (d-vectors)
    ! under symmetry operations of the parent multilattice
    ! are used later on. Generate them here
    tempD = d ! Make a temporary copy of the d-set members
    call matrix_inverse(parLV,inv_parLV,err); if(err) stop "Inverse failed for d-mapping"

    do iD = 1, nD
       call bring_into_cell(d(:,iD),inv_parLV,parLV,eps)
       if(.not. equal(d(:,iD),tempD(:,iD),eps)) then
          write(*,'("Original d-set member was not inside the unit cell. It has been remapped.")')
          write(*,'("Original:",3(f7.3,1x))') tempD(:,iD)
          write(*,'("Remapped:",3(f7.3,1x))') d(:,iD)
       endif
    enddo
    call get_dvector_permutations(parLV,d,ParRPList,LatDim,eps)
    
    ! This part generates all the derivative structures. Results are
    ! written to unit 14.
    Tcnt = 0 ! Keep track of the total number of structures generated
    HNFcnt = 0 ! Keep track of the total number of
    ! symmetrically-inequivalent HNFs in the output

    do ivol = nMin, nMax
       if (conc_check) then
          call get_list_at_this_size(ivol,nD,cRange,iRange,eps)
          if (size(iRange,1)==0) cycle
       endif
       call cpu_time(tstart)
       if (fixed_cells) then
          call read_in_cells_from_file(ivol,HNF,parLV,eps)
          if(size(HNF,3)==0) cycle
       else
          if (LatDim==3) then !<<< 2D ? or 3D? >>
             call get_all_HNFs(ivol,HNF)    ! 3D
          else
             call get_all_2D_HNFs(ivol,HNF) ! 2D
          endif
       endif
       
       !call cpu_time(HNFtime) Many of the superlattices will be
       ! symmetrically equivalent so we use the symmetry of the parent
       ! multilattice to reduce the list to those that are
       ! symmetrically distinct.
       call remove_duplicate_lattices(HNF,LatDim,parLV,d,ParRPList,rdHNF,fixOp,RPList,uqlatts,hnf_degen,eps)
       !call cpu_time(Removetime) Superlattices with the same SNF will
       ! have the same list of translation permutations of the
       ! labelings. So they can all be done at once if we find the
       ! SNF. rdHNF is the reduced list.
       call get_SNF(rdHNF,L,SNF,R,RPList,uqSNF,SNF_labels,fixOp)
       ! call cpu_time(SNFtime) Each HNF will have a certain number of
       ! rotations that leave the superlattice fixed, called
       ! fixOp. These operations will effect a permutation on the
       ! (d,g) table. Since many of the HNFs will have an identical
       ! list of rotation permutations, it'll be efficient to reduce
       ! the labelings for all such HNFs just once. So we need to
       ! generate the list for each HNF and then sort the HNFs into
       ! blocks with matching permutation lists. And we'll go ahead
       ! and include the translation permutations as well so that each
       ! list in rdRPList contains all possible permutations that
       ! identify duplicate labelings.

       call get_rotation_perms_lists(parLV,rdHNF,L,SNF,fixOp,RPList,ParRPList,eps,aperms,use_arrows=arrows,surf=(Latdim==2))
       call organize_rotperm_lists(RPList,rdRPList,RPLindx,aperms,rdaperms)
       
       ! This next if statement makes the run-time horrible (N^3
       ! scaling) if enabled. (only used for checking once.)
       Scnt = 0 ! Keep track of the number of structures at this size   
       do iBlock = 1, maxval(RPLindx)
          !call cpu_time(blockstart)
          ! filename = "debug_temp_perms.out"

          ! call write_rotperms_list(rdRPList(iBlock),filename)
          if (conc_check) then
             do iC = 1, size(iRange,1) ! loop over each concentration in the range
                ! If we are running in polya mode then we only want to
                ! run the polya code. Otherwise run the full enumeration.
                if (present(polya) .and. (polya .eqv. .true.)) then
                   call polya_count(iRange(iC,:), rdRPList(iBlock)%perm, rdaperms(iBlock)%perm,&
                        Tcnt, Scnt, iBlock, RPLindx)
                else
                   allocate(site_res(size(rdRPList(iBlock)%perm(1,:)),size(iRange(iC,:))))
                   site_res = 0
                   do i = 1, size(iRange(iC,:)); do j = 1, size(rdRPList(iBlock)%perm(1,:))
                      if (any(labelFull(:,(j-1)/ivol+1)==i-1)) then
                         site_res(j,i) = 1
                      end if
                   end do; end do

                   ! If there are site restrictions present in the system
                   ! we want to use the old enum3 code since it is more
                   ! efficient. Unless there are arrows present in the
                   ! enumeration or the multinomial of possible
                   ! arrangements is to large for enum3 to handle.
                   if (any(site_res == 0) .and. (multinomial(iRange(iC,:)) < max_binomial) .and. (arrows .eqv. .false.)) then
                      call generate_permutation_labelings(k,ivol,nD,rdRPList(iBlock)%perm,&
                           lm,iRange(iC,:),labelFull,digitFull,lab_degen,fixed_cells)
                      call write_labelings(k,ivol,nD,label,digit,iBlock,rdHNF,SNF,L,fixOp,Tcnt,&
                           Scnt,HNFcnt,RPLindx,lm,equivalencies,hnf_degen,lab_degen,iRange(iC,:))
                      
                   else
                      call recursively_stabilized_enum(rdRPList(iBlock)%perm,iRange(iC,:),ivol,k,SNF,L,rdHNF,HNFcnt,&
                        hnf_degen,Tcnt,Scnt,fixOp,iBlock,equivalencies,RPLindx,site_res,&
                        fixed_cells,rdaperms(iBlock)%perm)
                   end if
                   deallocate(site_res)
                end if
                ! call generate_disjoint_permutation_labelings(k,ivol,nD&
                !     &,rdRPList(iBlock)%perm,lm,iRange(iC,:),labelFull,digitFull,2)
             enddo
          else
             call generate_unique_labelings(k,ivol,nD,rdRPList(iBlock)%perm,&
                  full,lm,label,digit,lab_degen,fixed_cells)
             ! Now that we have the labeling marker, we can write the output.
             call write_labelings(k,ivol,nD,label,digit,iBlock,rdHNF,SNF,L,fixOp,Tcnt,Scnt,HNFcnt,RPLindx,lm,equivalencies,hnf_degen,lab_degen)
             !call cpu_time(endwrite)
          endif
       enddo! iBlock
       !call cpu_time(writetime)
       call cpu_time(tend)
       write(*,'(i4,1x,f14.4,1x,i8,3x,i3,3x,i7,7x,f7.4,i12,i12)')ivol,tend-tstart,size(HNF,3),&
            size(uqSNF,3),size(rdHNF,3),1-size(rdHNF,3)/real(size(HNF,3)),Scnt, Tcnt
    enddo ! loop over cell sizes (ivol)
    close(14)
    write(*,'(A)') "---------------------------------------------------------------------------------------------"
                
  ENDSUBROUTINE gen_multilattice_derivatives

  !!<summary>This subroutine finds the polya prediction for the number
  !!of unique arrangements at each unique concentration range and cell
  !!size.</summary>
  !!<parameter name="total_count" regular="true">The total number of
  !!structures predicted so far.</parameter>
  !!<parameter name="size_count" regular="true">The total number of
  !!structures predicted for this cell size.</parameter>
  !!<parameter name="concs" regular="true">The concentration we want
  !!the number of configurations for.</parameter>
  !!<parameter name="siteperms" regular="true">The permutation group
  !!for the sites.</parameter>
  !!<parameter name="arrowperms" regular="true">The permutations of
  !!the arrows.</parameter>
  SUBROUTINE polya_count(concs, siteperms, arrowperms, total_count, size_count, HNFi, permIndx)
    integer, pointer, intent(in) :: siteperms(:,:), arrowperms(:,:)
    integer, intent(in) :: concs(:), permIndx(:)
    integer, intent(inout) :: total_count, size_count
    integer, intent(in) :: HNFi

    !!<local name="poly">An array used by the polya algorithm to keep
    !!track of the polynomials.</local>
    !!<local name="status">Allocation status flag.</local>
    !!<local name="tconc">A copy of the concentrations that allow them
    !!to be sorted.</local>
    !!<local name="conc_map">Stores the mapping from the arrow labels
    !!to the non arrow labels</local>
    !!<local name="use_arrows">Logical, true if arrows.in is present.</local>
    !!<local name="this_count">The number of unique arrangnments
    !!predicted for this system.</local>
    !!<local name="labels">The labels of the atoms present.</local>
    !!<local name="species_i">Variable for loops.</local>
    !!<local name="species_j">Variable for loops.</local>
    !!<local name="a_conc">The concentrations with the arrows included.</local>
    !!<local name="arrows">The number of arrows for each color.</local>
    !!<local name="nArrows">The number of sites with arrows on them.</local>
    integer, allocatable :: poly(:,:), conc_map(:,:), poly_a(:,:,:)
    logical :: use_arrows
    integer :: status, species_i, species_j, nHNF, nArrows
    integer(li) :: this_count
    integer, allocatable :: tconc(:), labels(:), a_conc(:)
    integer :: arrows(size(concs))

    inquire(FILE="arrows.in",EXIST=use_arrows)
    ! Get the arrow concentrations and adjust the input concentrations
    ! if needed. Also create a mapping that will later allow us to
    ! undo the transformations in the colorings that happen in the
    ! next step.
    if (use_arrows) then
       call read_arrows(size(concs), arrows)
       call arrow_concs(concs,arrows,a_conc,conc_map,nArrows)
       allocate(poly(size(siteperms,1),size(siteperms,2)), poly_a(size(siteperms,1),size(siteperms,2),2), STAT=status)
       if(status/=0) stop "Allocation failed in recursively_stabilized_enum: poly poly_a."
    else
       allocate(a_conc(size(concs)), poly(size(siteperms,1),size(siteperms,2)), poly_a(size(siteperms,1),size(siteperms,2),2), STAT=status)
       if(status/=0) stop "Allocation failed in recursively_stabilized_enum: a_conc, poly poly_a."
       a_conc = concs
       allocate(conc_map(1,2),STAT=status)
       if(status/=0) stop "Allocation failed in recursively_stabilized_enum: conc_map."
       conc_map(1,:) = (/0,0/)
    end if
    allocate(tconc(count(a_conc > 0)),labels(count(a_conc > 0)),STAT=status)
    if(status/=0) stop "Allocation failed in recursively_stabilized_enum: tconc, poly, labels."

    ! remove any of the zero concentration elements from the list.
    species_j = 1
    do species_i = 1, size(a_conc)
       if (a_conc(species_i) > 0) then
          tconc(species_j) = a_conc(species_i)
          labels(species_j) = species_i-1
          species_j = species_j + 1
       end if
    end do

    ! sort the concentrations to be in the optimal order
    call heapsort(tconc,labels,conc_map) 

    ! Now we can finally run the polya algorithm.
    if ((size(tconc) == 1) .and. (.not. use_arrows)) then
       this_count = 1
    else if (.not. use_arrows) then
       this_count = polya(tconc, siteperms, polynomials=poly, decompose=.True.)
    else
       this_count = polya(tconc, siteperms, agroup=arrowperms, polynomials=poly_a, arrows=count(conc_map(:,1) >0), decompose=.True.)
    end if

    nHNF = count(permIndx == HNFi)
    size_count = size_count + this_count*nHNF
    total_count = total_count + this_count*nHNF
    deallocate(tconc,labels,a_conc,conc_map)
    if (allocated(poly)) deallocate(poly)
    if (allocated(poly_a)) deallocate(poly_a)
  end SUBROUTINE polya_count
  
END MODULE derivative_structure_generator
