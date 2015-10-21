!!<summary>Defines a trees class to be used for enumerating derivative structures at a fixed
!!concentration.</summary>
MODULE tree_class
  use combinatorics
  use group_theory
  use enumeration_types
  implicit none
  private
  public tree, enumerate_unique_permutations


  type :: GroupList
     type(permList), pointer :: layer(:) => null()
  endtype GroupList

  !!<summary>A "tree" object that contains all the members needed for
  !!enumerating derivative structures at a fixed
  !!concentration. Translated from the python version "enum4" in the
  !!projects git repository</summary>
  type :: tree 
     private
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
     !!<member name="done"> The entire tree has been traversed
     !!</member>
     integer, pointer :: colors(:) => null()
     integer          :: k 
     integer          :: n 
     type(GroupList)  :: G
     integer, pointer :: Gsize(:) => null()
     integer, pointer :: loc(:) => null()
     integer, pointer :: branches(:) => null()
     logical          :: done = .false.
   contains
     procedure, public :: init => initializeTree
     procedure, public :: coloring => generateColoringFromLocation 
     procedure, public :: ccIdx => indexOfCurrentColor
     procedure, public :: depth 
     procedure, public :: increment_location
     procedure, public :: next_branch
     procedure, public :: enumerate_unique_permutations
     procedure, public :: listGroup
  endtype tree
  

  !!<summary>A list of group operations for each layer in the
  !!tree. I.e., the list of stabilizer groups</summary>

  !!<summary>Stores a list of permutations (for the atomic sites in
  !!deriv. structure)</summary>
  type :: permList
     integer, pointer :: perms(:,:) => null()
  endtype permList

  
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
  !!<parameter name="colors">A list of the numbers of each color in
  !!the enumeration</parameter>
  !!<parameter name="generators">Symmetry operators (permutations)
  !!that generate the group (but may contain the whole
  !!group)</parameter>
  !!<parameter name="makeG">Flag: True-> Generate group, False ->
  !!Generators already make a group</parameter>
  subroutine initializeTree(self,colors,generators,makeG)
    class(tree)         :: self
    integer, intent(in) :: colors(:)
    integer, pointer    :: generators(:,:)
    logical, optional   :: makeG
    
    integer i

    self%k = size(colors) 
    allocate(self%colors(self%k))
    self%colors = colors
    self%n = sum(colors) 
    allocate(self%loc(self%k),self%branches(self%k-1))
    allocate(self%Gsize(self%k))
    self%loc = -1
    do i = 1, self%k-1
       self%branches(i) = nchoosek(sum(self%colors(i:)),self%colors(i))
    enddo
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
       else 
          allocate(self%G%layer(1)%perms(size(generators,1),self%n))
       endif
    else
       allocate(self%G%layer(1)%perms(size(generators,1),self%n))
    endif
    self%Gsize(3:self%k) = 0
    self%Gsize(1:2) = size(self%G%layer(1)%perms,1)
  endsubroutine initializeTree

  !!<summary>Generate the coloring associated with the current location
  !!in the tree</summary>
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
  
  !!<summary>Return the index for the current color (i.e., the current
  !!layer in the tree)</summary>
  function indexOfCurrentColor(self)
    class(tree), intent(in) :: self
    integer    :: indexOfCurrentColor
    
    indexOfCurrentColor = minloc(self%loc,1,self%loc==-1) -1 ! Current
    ! location has -1's for all layers below current depth
  end function indexOfCurrentColor


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
    integer             :: integer2coloring(m)
    integer, intent(in) :: y,m,a

    integer I,t,ell
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
    ! depth is 1 less that location of first -1
  endfunction depth

  !!<summary>Increment the location in tree. Either move across
  !!branches at the same depth or move up or down between levels. When
  !!possible downward moves are first, then lateral moves, lastly
  !!moves up the tree.</summary>
  !!<parameter name="self" regular="true"></parameter>
  subroutine increment_location(self)
    class(tree) :: self

    integer :: d,i

    d = self%depth()
    if (d < self%k - 1) then ! we can still go down in the tree
       d = d + 1
       self%loc(d) = 0
       allocate(self%G%layer(d+1)%perms(self%Gsize(d),self%n)) ! Set
       ! next stabilizer to be as big as previous
    else ! We are at the bottom of the tree (2nd lowest layer, but
       ! lowest always has just one coloring)
       self%loc(d) = self%loc(d) + 1
!GLWH> Don't need to reset this stabilizer because we are on the lowest layer?
!   self%Gsize(d+1) = 0 ! Reset the number of stabilizer elements for next level down
       do while (self%loc(d) >= self%branches(d))! If at the end of branches on this layer, move
          d = d - 1                            ! up until we can go down again
          if (d < 1) then
             self%done = .true.
             exit ! All done with the tree
          endif
          deallocate(self%G%layer(d+1)%perms) ! reset the stabilizer for this depth
          self%Gsize(d+1) = 0
          self%loc(d+1) = -1
          self%loc(d) = self%loc(d) + 1
       enddo
       if (d > 0) allocate(self%G%layer(d+1)%perms(self%Gsize(d),self%n))
    endif
  endsubroutine increment_location

  !!<summary>Advance the tree to the next branch. Move up to a higher
  !!level if the branches at this level are exhausted</summary>
  !!<parameter name="self" regular="true"></parameter>
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
       ! elements at this level?  Use that at next (upper limit)
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
       self%loc(d) = self%loc(d) + 1
       self%Gsize(d+1) = 0 ! Reset the counter for the stabilizer ops
       allocate(self%G%layer(d+1)%perms(self%Gsize(d),self%n))
    enddo
  endsubroutine next_branch
  
  
  subroutine listGroup(self)
    class(tree) :: self
    integer :: i
    
    write(*,'("Group size: ",2(i4,1x))') shape(self%G%layer(1)%perms)
    do i = 1, size(self%G%layer(1)%perms,1)
       write(*,'(i3,1x,100(i2,1x))') i, self%G%layer(1)%perms(i,:)
    enddo
    
  endsubroutine listGroup
END MODULE tree_class
