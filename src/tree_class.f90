!!<summary>Defines a trees class to be used for enumerating derivative structures at a fixed
!!concentration.</summary>
MODULE tree_class
  use combinatorics
  use group_theory
  use enumeration_types
  implicit none
  private
  public tree, enumerate_unique_permutations


  !!<summary>A "tree" object that contains all the members needed for
  !!enumerating derivative structures at a fixed
  !!concentration. Translated from the python version "enum4" in the
  !!projects git repository</summary>
  type :: tree 
     !!<member name="done"> The entire tree has been traversed
     !!</member>
     !!<member name="unique">A logical that indicates if a
     !!configuration is unique.</member>
     !!<member name="colors"> Number of each type (color); i.e.,
     !!number of layers in the tree </member>
     !!<member name="k"> Number of different colors </member>
     !!<member name="n"> Total number of objects (number of "slots")
     !!in the enumeration </member>
     !!<member name="G"> Symmetry groups for each layer in the tree
     !!</member>
     !!<member name="Gsize"> Counter for the number of operations in
     !!each stabilizer (group at each layer)</member>
     !!<member name="loc"> Location in the tree. k components, each
     !!component the branch number</member>
     !!<member name="branches"> Number of branches at each layer in
     !!the tree (k-1 components; last layer not needed) </member>
     !!<member name="base">An array that keeps track of the unique
     !!arrangements of the first layer of the enumeration.</member>
     integer, pointer :: colors(:) => null()
     integer          :: k 
     integer          :: n 
     type(GroupList)  :: G
     integer, pointer :: Gsize(:) => null()
     integer, pointer :: loc(:) => null()
     integer, pointer :: branches(:) => null()
     integer, pointer :: base(:) => null()
     logical          :: unique = .true.
     logical          :: done = .false.
   contains
     procedure, public :: init => initializeTree
     procedure, public :: coloring => generateColoringFromLocation 
     procedure, public :: depth 
     procedure, public :: increment_location
     procedure, public :: next_branch
     procedure, public :: enumerate_unique_permutations
     procedure, public :: listGroup
     procedure, public :: check => check_labeling
     procedure, public :: get_loc => generateLocationFromColoring
  endtype tree
  

CONTAINS
  !!<summary>Enumerate symmetrically-distinct labelings using the
  !!"recursive stabilizer" approach (enum4).</summary>
  !!<parameter name="self" regular="true"></parameter>
  !!<parameter name="colors" regular="true">A list of the numbers of
  !!each color in the enumeration</parameter>
  !!<parameter name="generators" regular="true">Symmetry operators
  !!(permutations) that generate the group (but may contain the whole
  !!group)</parameter>
  !!<parameter name="makeG" regular="true">Flag: True-> Generate
  !!group, False -> Generators already make a group</parameter>
  !!<parameter name="uqperms" regular="true">Output. List of the
  !!symmetrically-distinct permutations.</parameter>
  subroutine enumerate_unique_permutations(self,colors,generators,makeG,uqperms)
    class(tree)         :: self
    integer, intent(in) :: colors(:)
    integer, pointer    :: generators(:,:) 
    logical             :: makeG
    integer, pointer    :: uqperms(:,:) ! Output from the routine,
    ! list of the permutations that are distinct
    
    integer, allocatable:: labeling(:), rlabeling(:) ! Copy of
    ! labeling for current location and rotated (permuted) version of
    ! the same labeling according to one of the symmetry operations
    integer, allocatable:: tempPerms1(:,:), tempPerms2(:,:)

    integer             :: ig, cg ! Loop counter for the group elements, next stabilizer counter
    integer             :: i ! generic counter
    logical             :: dup ! Flag to indicate if a labeling is a duplicate
    integer             :: cuq ! Number of unique permutations (survivors) found so far

    call self%init(colors,generators,.true.)
    allocate(tempPerms1(100,self%k)) ! Start with 100 possible
    ! survivors (doubles each time this isn't enough)
    allocate(labeling(self%n),rlabeling(self%n))
    call self%increment_location()
    cuq = 0
    !do while (any(self%loc(self%k-1)<self%branches))
    do while (.not. self%done)
       call self%coloring(labeling) ! Copy the current labeling (the
       ! coloring for current location in the tree) Apply the group
       ! operations to this labeling. If the index of the rotated
       ! labeling is less than the index of the unrotated labeling,
       ! then it has already been hit in the list. It's a
       ! symmetrically-equivalent duplicate shouldn't be included in
       ! the list.
       dup = .false.
       !   write(*,'("<< loc >> ",10(i5,1x))') self%loc
       !   write(*,'(" || depth || ",i5)') self%depth()
       self%Gsize(self%depth()+1) = 0
       groupCheck: do ig = 1, self%Gsize(self%depth())
          ! Use the ig-th permutation to rearrarange "labeling" (i.e.,
          ! the permutation is a vector subscript)
          rlabeling = labeling(self%G%layer(self%depth())%perms(ig,:))
          if (all(rlabeling == labeling)) then ! this group op fixes the
             ! current labeling and is part of the stabilizer group
             self%Gsize(self%depth()+1) = self%Gsize(self%depth()+1) + 1
             cg = self%Gsize(self%depth()+1) 
             self%G%layer(self%depth()+1)%perms(cg,:) = self%G%layer(self%depth())%perms(ig,:)
          else ! Check to see if it is a duplicate
             do i = 1, self%n
                if (labeling(i) < rlabeling(i)) then ! this rlabeling
                ! appeared earlier in the list
                   call self%next_branch()
                   dup = .true.
                   exit groupCheck
                elseif (labeling(i) > rlabeling(i)) then ! this rlabeling is later in the list
                   exit
                endif
             enddo
          endif
       enddo groupCheck
       if (.not. dup) then ! if the labeling is a complete one, it is a survivor (unique)
          if (self%depth() == self%k - 1) then ! Labeling is complete
             cuq = cuq + 1
             if (cuq > size(tempPerms1,1)) then ! Need more storage, so expand the survivor list
                allocate(tempPerms2(cuq*2,self%k))
                tempPerms2(:size(tempPerms1,1),:self%k) = tempPerms1
                call move_alloc(tempPerms2,tempPerms1) ! Copy and deallocate temp variable
             endif
             tempPerms1(cuq,:) = self%loc 
          endif
          call self%increment_location()
       endif
    enddo
    ! Allocate uqperms to the actual number of survivors
    allocate(uqperms(cuq,self%k-1))
    uqperms = tempPerms1(1:cuq,:self%k-1)
  endsubroutine enumerate_unique_permutations

  !!<summary>Initialize the tree using "colors" for input</summary>
  !!<parameter name="colors" regular="true">A list of the numbers of each color in
  !!the enumeration</parameter>
  !!<parameter name="generators" regular="true">Symmetry operators (permutations)
  !!that generate the group (but may contain the whole
  !!group)</parameter>
  !!<parameter name="makeG" regular="true">Flag: True-> Generate group, False ->
  !!Generators already make a group</parameter>
  subroutine initializeTree(self,colors,generators,makeG)
    class(tree)         :: self
    integer, intent(in) :: colors(:)
    integer, pointer, intent(in)    :: generators(:,:)
    logical, optional, intent(in)   :: makeG
    
    integer :: i

    self%k = size(colors) 
    allocate(self%colors(self%k))
    self%colors = colors
    self%n = sum(colors) 
    allocate(self%loc(self%k),self%branches(self%k-1))
    allocate(self%Gsize(self%k))
    self%loc = -1
    if (self%k > 1) then
       allocate(self%branches(self%k-1))
       do i = 1, self%k-1
          self%branches(i) = nchoosek(sum(self%colors(i:)),self%colors(i))
       enddo
       allocate(self%base(self%branches(1)))
       self%base = 0
       ! Initialize the group for layer 1. Generate the group if makeG =
       ! .true., otherwise assume that the list of generators given is
       ! the entire group
       allocate(self%G%layer(self%k)) !Allocate number of perm
       !lists. Don't need one for last color but makes the code
       !cleaner. Perhaps rethink this later
       allocate(self%G%layer(1)%perms(size(generators,1),size(generators,2)))
       self%G%layer(1)%perms = generators
    
       if (present(makeG)) then
          if (makeG) then
             call grouper(self%G%layer(1)%perms)
          end if
       end if
       self%Gsize(3:self%k) = 0
       self%Gsize(1:2) = size(self%G%layer(1)%perms,1)
    else if (self%k < 1) then
       stop "ERROR in initialize tree. There are less than one atomic species present."
    end if

  endsubroutine initializeTree

  !!<summary>Generate the coloring associated with the current location
  !!in the tree</summary>
  !!<parameter name="labeling" regular="true">The output labeling for
  !!the location.</parameter>
  subroutine generateColoringFromLocation(self,labeling)
    class(tree) :: self
    integer, allocatable, intent(out) :: labeling(:)
    
    integer, allocatable :: clabeling(:), freeIndices(:), configList
    integer              :: ik, cIdx, iIdx, jIdx, ilc, jlc, nEmp
    
    allocate(labeling(self%n)); labeling = 0 ! Empty labeling container to start
    do ik = 1, self%depth() ! Loop over all colors up to the current depth
       nEmp = count(labeling==0) ! How many empty slots remaining up to this level?
       allocate(freeIndices(nEmp))
       jlc = 0 ! Counter for the empty slots
       ! Find the indices of the empty slots
       do ilc = 1, self%n
          if (labeling(ilc) == 0) then
             jlc = jlc + 1
             freeIndices(jlc) = ilc
          endif
       enddo
       cIdx = self%loc(ik) ! Index of the current color
       allocate(clabeling(nEmp))
       ! Get the partial coloring for the current color index so it can be added to the labeling
       clabeling = integer2coloring(cIdx,nEmp,self%colors(ik))
       ! Load up the empty slots that need the current colors
       do iIdx = 1, nEmp
          if (clabeling(iIdx)/=0) then
             labeling(freeIndices(iIdx)) = ik
          endif
       enddo
       deallocate(freeIndices,clabeling)
    enddo
  end subroutine generateColoringFromLocation
  
  !!<summary>Convert an integer to binary labeling. Takes the index of
  !!the current coloring, number of open slots, and number of this
  !!color. There are n (m) choose k (a) possible configurations and
  !!the integer (y) is the index in that list of
  !!configurations.</summary>
  !!<comments>Follows the algorithm in the enum3 paper, Comp Mat Sci
  !!59 101 (2010) exactly</comments>
  !!<parameter name="y" regular="true"></parameter>
  !!<parameter name="m" regular="true"></parameter>
  !!<parameter name="a" regular="true"></parameter>
  function integer2coloring(y,m,a)
    integer, intent(in) :: y,m,a
    integer             :: integer2coloring(m)

    integer :: I,t,ell
    integer, allocatable :: configList(:)
    
    I = y; t = a; ell = m
    allocate(configList(0:m-1)); configList = -1
    do while (ell > 0)
       if (nchoosek(ell-1,t-1) <= I) then
          configList(m-ell) = 0
          I = I - nchoosek(ell-1,t-1)
       else
          configList(m-ell) = 1
          t = t - 1
       endif
       ell = ell - 1
    enddo
    if (any(configList==-1)) stop "Error in integer2coloring: -1's remaining in configList"
    integer2coloring = configList
  endfunction integer2coloring

  !!<summary>Return the depth of the current location in the
  !!tree. Depth is indexed from 1 to k.</summary>
  !!<parameter name="self" regular="true"></parameter>
  function depth(self)
    integer     :: depth
    class(tree) :: self
    depth = minloc(self%loc,1,self%loc==-1)-1 ! Subtract 1 because the
    if (depth == -1) then
       depth = size(self%loc,1)
    end if
    ! depth is 1 less that location of first -1
  endfunction depth
  
  !!<summary>Increment the location in tree. Either move across
  !!branches at the same depth or move up or down between levels. When
  !!possible downward moves are first, then lateral moves, lastly
  !!moves up the tree.</summary>
  subroutine increment_location(self)
    class(tree) :: self

    !!<local name="skipping">Indicates if we're still skipping
    !!configurations on the base level.</local>
    integer :: d, i
    logical :: skipping

    d = self%depth()
    if ((d < self%k - 1) .and. (self%unique .eqv. .true.)) then ! we can still go down in the tree
       d = d + 1
       self%loc(d) = 0
       
       ! Set next stabilizer to be as big as previous
       allocate(self%G%layer(d+1)%perms(self%Gsize(d),self%n)) 
       self%G%layer(d+1)%perms = 0
    else ! We are at the bottom of the tree (2nd lowest layer, but
       ! lowest always has just one coloring)
       self%loc(d) = self%loc(d) + 1
       if (self%loc(d) >= self%branches(d)) then
          deallocate(self%G%layer(d+1)%perms) ! reset the stabilizer for the next layer
          self%Gsize(d+1) = 0
          do while (self%loc(d) >= self%branches(d))! If at the end of branches on this layer,
                                                   ! move up until we can go down again
             d = d - 1
             if (d < 1) then
                self%done = .true.
                exit ! All done with the tree
             endif
             deallocate(self%G%layer(d+1)%perms) ! reset the stabilizer for this depth
             self%Gsize(d+1) = 0
             self%loc(d+1) = -1
             self%loc(d) = self%loc(d) + 1
          end do
       end if
      
       if (d == 1) then
          skipping = .true.
          do while (skipping .eqv. .true.)
             if (self%loc(d) >= self%branches(d)) then
                self%done = .true.
                skipping = .false.
                d = d -1
             else if (self%base(self%loc(d)+1) == 0) then
                skipping = .false.
             else
                self%loc(d) = self%loc(d) + 1
             end if
          end do
       end if

       if (d > 0) then
          allocate(self%G%layer(d+1)%perms(self%Gsize(d),self%n))
          self%G%layer(d+1)%perms = 0
       end if
    endif
  endsubroutine increment_location

  !!<summary>Advance the tree to the next branch. Move up to a higher
  !!level if the branches at this level are exhausted</summary>
  subroutine next_branch(self)
    class(tree) :: self
    integer :: d, nSg

    d = self%depth()
    self%loc(d) = self%loc(d) + 1 ! Move one branch to the right
    self%Gsize(d+1) = 0 ! Reset the counter for the number of group ops in the stabilizer
    !GLWH Not sure the next if block really does anything. Rethink this when the rest is working
    if (associated(self%G%layer(d+1)%perms)) then
       deallocate(self%G%layer(d+1)%perms) ! Reset the stabilizer one level below here
       nSg = size(self%G%layer(d)%perms,1) ! How many stabilizer
       ! elements at this level?  Use that as next (upper limit)
       allocate(self%G%layer(d+1)%perms(nSg,self%n))
    endif
    do while (self%loc(d) >= self%branches(d)) ! We are at the end of
    ! the branches at this level so go up
       d = d - 1
       self%loc(d+1) = -1
       if (d < 1) then ! All done with the tree, exit loop
          self%done = .true.
          exit 
       endif
       deallocate(self%G%layer(d+1)%perms)
       if (d == 1) then
          do while (self%base(self%loc(d)) == 1)
             self%loc(d) = self%loc(d) + 1
          end do
       else
          self%loc(d) = self%loc(d) + 1
       end if
       self%Gsize(d+1) = 0 ! Reset the counter for the stabilizer ops
       allocate(self%G%layer(d+1)%perms(self%Gsize(d),self%n))
    enddo
  endsubroutine next_branch
  
  !!<summary>Writes the group to the screen.</summary>
  subroutine listGroup(self)
    class(tree) :: self
    integer :: i
    
    write(*,'("Group size: ",2(i4,1x))') shape(self%G%layer(1)%perms)
    do i = 1, size(self%G%layer(1)%perms,1)
       write(*,'(i3,1x,100(i2,1x))') i, self%G%layer(1)%perms(i,:)
    enddo
    
  endsubroutine listGroup

  !!<summary>This routine checks to see if a labeling is unique given
  !!a permutation group.</summary>
  !!<parameter name="label" regular="true">The labeling that is
  !!being checked.</parameter>
  !!<parameter name="symsize" regular="true">Then number of cells being used
  !!in the system.</parameter>
  !!<parameter name="fixedcell" regular="true">A logical for if this is a
  !!fixed cell.</parameter>
  SUBROUTINE check_labeling(self,label,symsize,fixedcell)
    class(tree)         :: self
    integer, intent(in) :: label(:)
    integer, intent(in) :: symsize
    logical, intent(in) :: fixedcell

    !!<local name="rotatedlabel">The labeling after a permutation
    !!has been applied.</local>
    !!<local name="i">Indexing variable</local>
    !!<local name="j">Indexing variable</local>
    !!<local name="stab_size">The size of the next stabilizer group</local>
    !!<local name="new_loc">The location of the rotated branch.</local>
    !!<local name="min_stab">Temporary array for the stabilizers</local>
    integer :: rotatedlabel(size(label))
    integer :: i, j, stab_size, new_loc
    integer, allocatable :: min_stab(:,:)
    
    ! Make sure that the number of stibalizers is initially 0
    self%Gsize(self%depth() + 1) = 0

    ! Now apply the symmetry group in order to determine if the
    ! labeling is unique and to find the stabilizers
    groupCheck: do i = 1, self%Gsize(self%depth())
       ! Apply the ith permutation to the labeling
       rotatedlabel = label(self%G%layer(self%depth())%perms(i,:))
       if (all( rotatedlabel == label )) then
          if ((self%depth() == self%k -1).and. (i <= symsize) .and. (i > 1) .and. (.not. fixedcell)) then
             do j = 2, symsize
                if (all(self%G%layer(1)%perms(j,:) == self%G%layer(self%depth())%perms(i,:))) then
                   self%unique = .False.
                   exit groupCheck
                end if
             end do
          end if
          ! This permutation is a stabilizer, i.e., it fixes the
          ! current labeling
          stab_size = self%Gsize(self%depth()+1) + 1
          self%Gsize(self%depth()+1) = stab_size
          self%G%layer(self%depth()+1)%perms(stab_size,:) = self%G%layer(self%depth())%perms(i,:)
       else 
          ! check to see if the configuration is unique
          call self%get_loc(rotatedlabel,new_loc)
          if (self%depth() == 1) then
             self%base(self%loc(self%depth())+1) = 1
             if (new_loc < self%loc(self%depth())) then
                ! This labeling appeared earlier in the list
                self%unique = .False.
                exit groupCheck
             elseif (new_loc > self%loc(self%depth())) then
                self%base(new_loc + 1) = 1
             end if
          else
             if (new_loc < self%loc(self%depth())) then
                ! The rotated labeling is equivalent
                ! to an earlier configuration.
                self%unique = .False.
                exit groupCheck
             end if
          end if
       end if
    end do groupCheck

    if (self%unique .eqv. .True.) then
       allocate(min_stab(stab_size,self%n))
       min_stab = self%G%layer(self%depth()+1)%perms(1:stab_size,:)
       deallocate(self%G%layer(self%depth()+1)%perms)
       allocate(self%G%layer(self%depth()+1)%perms(stab_size,self%n))
       self%G%layer(self%depth()+1)%perms = min_stab
       deallocate(min_stab)
    end if
   
  END SUBROUTINE check_labeling

  !!<summary>This routine takes a labeling and returns it's location
  !!in the branch of the tree.</summary>
  !!<parameter name="labeling" regular="true">The labeling whose
  !!location is to be determined.</parameter>
  !!<parameter name="index" regular="true">The output location.</parameter>
  SUBROUTINE generateLocationFromColoring(self,labeling,index)
    class(tree)          :: self
    integer, intent(in)  :: labeling(:)
    integer, intent(out) :: index

    !!<local name="new_labeling">Contains a binary labeling, i.e.,
    !!'1's and '0's.</local>
    !!<local name="i">Counting term index.</local>
    !!<local name="j">Counting term index.</local>
    !!<local name="nzeros">The number of zeros left in the labeling.</local>
    !!<local name="index">The output location.</local>
    !!<local name="ndigits">The number of digits to the current zero.</local>
    !!<local name="numones">The number of '1's to the right of the current zero.</local>
    !!<local name="minl">The location of the next zero.</local>
    !!<local name="lastone">The index of the last '1' in the binary array.</local>
    integer, allocatable :: new_labeling(:)
    integer :: i, j, nzeros, numones, lastone
    integer :: ndigits(1), minl(1)

    ! The labeling needs to have only '1's for the current color, have
    ! all smaller integer colorings removed, and be zeros
    ! otherwise. We also keey track of the location of the last '1' in the list.

    ! better algorithm for this if we exclude the looping over the concs
    !     do iK = 1, n
    !    allocate(mask(nLeft))
    !    mask = 0
    !    where(pack(l,l>=iK-1)==iK-1)
    !       !   where(l==iK-1)
    !       mask = 1
    !    end where
    !    xTemp = 0
    !    do iM = 1, nLeft
    !       if (mask(iM)==0) then!&
    !          xTemp = xTemp + binomial(nLeft - iM, count(mask(iM:)==1)-1) 
    !       endif
    !    enddo
    !    X(iK) = xTemp
    !    nLeft = nLeft - conc(iK)
    !    deallocate(mask)
    ! enddo

    allocate(new_labeling(count(labeling .eq. 0,1) + count(labeling .eq. self%depth(),1)))
    new_labeling = 0
    j = 1

    do i = 1, size(labeling)
       if (labeling(i) == self%depth()) then
          new_labeling(j) = 1
          lastone = j
          j = j + 1
       else if (labeling(i) == 0) then
          new_labeling(j) = 0
          j = j + 1
       end if
    end do

    ! Now we need the number of zeros between the beginning of the
    ! array and the last '1'.
    nzeros = 0!maxloc(new_labeling,1)
    do i = 1, lastone!, maxloc(new_labeling,1)
       if (new_labeling(i) == 0) then
          nzeros = nzeros + 1
       end if
    end do
    index = 0

    ! Use the binomial coefficient to find the index for this labeling.
    do i = 1, nzeros
       minl = minloc(new_labeling)
       ndigits = size(new_labeling) - minloc(new_labeling)
       numones = count(new_labeling(minl(1):) .eq. 1)
       new_labeling(minloc(new_labeling)) = 1
       index = index + binomial(ndigits(1),numones-1)
    end do

  END SUBROUTINE  generateLocationFromColoring
  
END MODULE tree_class
