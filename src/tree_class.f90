!!<summary>Defines a trees class to be used for enumerating derivative structures at a fixed
!!concentration.</summary>
MODULE tree_class
  use combinatorics
  use group_theory
  use enumeration_types
  use arrow_related
  implicit none
  private
  public tree, write_single_labeling


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
     !!<member name="nArrows">The integer number of arrows in the
     !!system.</member>
     !!<member name="color_map">The list of colors that have arrows
     !!associated with them.</member>
     !!<member name="A">The permutation group of the arrows.</member>
     integer, pointer :: colors(:) => null()
     integer          :: k 
     integer          :: n 
     type(GroupList)  :: G
     type(GroupList)  :: A
     integer, pointer :: Gsize(:) => null()
     integer, pointer :: loc(:) => null()
     integer, pointer :: branches(:) => null()
     integer, pointer :: base(:) => null()
     logical          :: unique = .true.
     logical          :: done = .false.
     integer          :: nArrows = 0
     integer, pointer :: color_map(:,:) => null()
   contains
     procedure, public :: init => initializeTree
     procedure, public :: coloring => generateColoringFromLocation 
     procedure, public :: depth 
     procedure, public :: increment_location
     procedure, public :: check => check_labeling
     procedure, public :: get_loc => generateLocationFromColoring
     procedure, public :: add_arrows => addArrowsToEnumeration
  endtype tree
  

CONTAINS

  !!<summary>Initialize the tree using "colors" for input</summary>
  !!<parameter name="colors" regular="true">A list of the numbers of each color in
  !!the enumeration</parameter>
  !!<parameter name="generators" regular="true">Symmetry operators (permutations)
  !!that generate the group (but may contain the whole
  !!group)</parameter>
  !!<parameter name="makeG" regular="true">Flag: True-> Generate site
  !!permutation group (will not generate arrow permutations), False ->
  !!Generators already make a group</parameter>
  !!<parameter name="arrow_group" regular="true">The arrow group for the system.</parameter>
  !!<parameter name="color_map" regular="true">The mapping of the arrowed colors back
  !! to the non-arrowed colors.</parameter>
  !!<parameter name="nArrows">The number of sites with arrows on them.</parameter>
  subroutine initializeTree(self,colors,generators, arrow_group, color_map, nArrows, makeG)
    class(tree)         :: self
    integer, intent(in) :: colors(:)
    integer, intent(in) :: color_map(:,:)
    integer, pointer, intent(in) :: arrow_group(:,:)
    integer, pointer, intent(in)    :: generators(:,:)
    integer, intent(in) :: nArrows
    logical, optional, intent(in)   :: makeG
    
    integer :: species_i, status

    self%k = size(colors)
    allocate(self%colors(self%k),STAT=status)
    if(status/=0) stop "Allocation failed in initializeTree: self%colors."       
    self%colors = colors
    self%n = sum(colors) 
    allocate(self%loc(self%k),self%branches(self%k-1),STAT=status)
    if(status/=0) stop "Allocation failed in initializeTree: self%loc, self%branches."       
    allocate(self%Gsize(self%k),STAT=status)
    if(status/=0) stop "Allocation failed in initializeTree: self%Gsize."           
    self%loc = -1
    self%nArrows = nArrows
    allocate(self%color_map(size(color_map,1),size(color_map,2)),STAT=status)
    if(status/=0) stop "Allocation failed in initializeTree: self%color_map."       
    self%color_map = color_map
    if (self%k > 1) then
       allocate(self%branches(self%k-1),STAT=status)
       if(status/=0) stop "Allocation failed in initializeTree: self%branches."       
       do species_i = 1, self%k-1
          self%branches(species_i) = nchoosek(sum(self%colors(species_i:)),self%colors(species_i))
       enddo
       allocate(self%base(self%branches(1)),STAT=status)
       if(status/=0) stop "Allocation failed in initializeTree: self%base"       
       self%base = 0
       ! Initialize the group for layer 1. Generate the group if makeG =
       ! .true., otherwise assume that the list of generators given is
       ! the entire group
       allocate(self%G%layer(self%k),self%A%layer(self%k),STAT=status) !Allocate number of perm lists.
       if(status/=0) stop "Allocation failed in initializeTree: self%G%layer, self%A%layer" 
       allocate(self%G%layer(1)%perms(size(generators,1),size(generators,2)),STAT=status)
       if(status/=0) stop "Allocation failed in initializeTree: self%G%layer%perms"       
       allocate(self%A%layer(1)%perms(size(arrow_group,1),size(arrow_group,2)),STAT=status)
       if(status/=0) stop "Allocation failed in initializeTree: self%A%layer%perms"       
       self%G%layer(1)%perms = generators
       self%A%layer(1)%perms = arrow_group

       ! makeG is a flag that tells the code it was passed generators
       ! for the group and not the entire group. This is useful for
       ! debugging the code so it has been kept in the release
       ! version. Please note that at present only the permutations of
       ! the sites can be constructed from the generators. This cannot
       ! be used to construct the arrow group.
       if (present(makeG)) then
          if (makeG) then
             call grouper(self%G%layer(1)%perms)
          end if
       end if
       self%Gsize(3:self%k) = 0
       self%Gsize(1:2) = size(self%G%layer(1)%perms,1)
    else if (self%k < 1) then
       stop "ERROR in initialize tree. There are less than one atomic species present."
    else if (self%k == 1) then
       ! self%done = .True.
       self%unique = .True.
       allocate(self%G%layer(self%k),self%A%layer(self%k),STAT=status) !Allocate number of perm lists
       if(status/=0) stop "Allocation failed in initializeTree: self%G%layer, self%A%layer" 
       allocate(self%G%layer(1)%perms(size(generators,1),size(generators,2)),STAT=status)
       if(status/=0) stop "Allocation failed in initializeTree: self%G%layer%perms"       
       allocate(self%A%layer(1)%perms(size(arrow_group,1),size(arrow_group,2)),STAT=status)
       if(status/=0) stop "Allocation failed in initializeTree: self%A%layer%perms"       
       self%G%layer(1)%perms = generators
       self%A%layer(1)%perms = arrow_group
    end if
  endsubroutine initializeTree

  !!<summary>Generate the coloring associated with the current location
  !!in the tree. This subroutine is an implementation of
  !!the algorithm described in section 3.1 of
  !!http://msg.byu.edu/papers/enum3.pdf (graphically displayed in Fig
  !!5).</summary>
  !!<parameter name="labeling" regular="true">The output labeling for
  !!the location.</parameter>
  subroutine generateColoringFromLocation(self,labeling)
    class(tree) :: self
    integer, allocatable, intent(out) :: labeling(:)
    
    integer, allocatable :: clabeling(:), freeIndices(:), configList
    integer              :: ik, cIdx, iIdx, jIdx, ilc, jlc, nEmp, status

    allocate(labeling(self%n),STAT=status)
    if(status/=0) stop "Allocation failed in generateColoringFromLocation: labeling."       
    labeling = 0 ! Empty labeling container to start
    do ik = 1, self%depth() ! Loop over all colors up to the current depth
       nEmp = count(labeling==0) ! How many empty slots remaining up to this level?
       allocate(freeIndices(nEmp),STAT=status)
       if(status/=0) stop "Allocation failed in generateColoringFromLocation: freeIndices." 
       jlc = 0 ! Counter for the empty slots
       ! Find the indices of the empty slots
       do ilc = 1, self%n
          if (labeling(ilc) == 0) then
             jlc = jlc + 1
             freeIndices(jlc) = ilc
          endif
       enddo
       cIdx = self%loc(ik) ! Index of the current color
       allocate(clabeling(nEmp),STAT=status)
       if(status/=0) stop "Allocation failed in generateColoringFromLocation: clabeling."
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

    integer :: I,t,ell, status
    integer, allocatable :: configList(:)
    
    I = y; t = a; ell = m
    allocate(configList(0:m-1),STAT=status)
    if(status/=0) stop "Allocation failed in integer2coloring: configList."       
    configList = -1
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
    !!<local name="d">Local storage of the current depth in the
    !!tree.</local>
    integer :: d, status
    logical :: skipping

    d = self%depth()

    if ((d < self%k - 1) .and. (self%unique .eqv. .true.)) then
       ! we can still go down in the tree
       d = d + 1
       self%loc(d) = 0
       
       ! Set next stabilizer to initially be as big as the previous layer
       allocate(self%G%layer(d+1)%perms(self%Gsize(d),self%n),STAT=status)
       if(status/=0) stop "Allocation failed in increment_location: self%G%layer%perms."
       allocate(self%A%layer(d+1)%perms(self%Gsize(d),size(self%A%layer(d)%perms,2)),STAT=status) 
       if(status/=0) stop "Allocation failed in increment_location: self%A%layer%perms."
       self%G%layer(d+1)%perms = 0
       self%A%layer(d+1)%perms = 0
    else if (self%k /= 1) then
       ! We are at the bottom of the tree (2nd lowest layer, but the
       ! lowest layer always has just one coloring) so we just need to
       ! move to the next branch at this layer and reset the
       ! stabilizers for the next layer (We keep track of the last
       ! layer's stabilizers so we can check for superperiodic
       ! structures in the check_labeling subroutine).
       self%loc(d) = self%loc(d) + 1
       if (self%loc(d) >= self%branches(d)) then
          deallocate(self%G%layer(d+1)%perms) ! reset the stabilizer for the next layer
          deallocate(self%A%layer(d+1)%perms) ! reset the stabilizer for the next layer
          self%Gsize(d+1) = 0
          do while (self%loc(d) >= self%branches(d))! If at the end of branches on this layer,
                                                   ! move up until we can go down again
             d = d - 1
             if (d < 1) then
                self%done = .true.
                exit ! All done with the tree
             endif
             deallocate(self%G%layer(d+1)%perms) ! reset the stabilizer for this depth
             deallocate(self%A%layer(d+1)%perms) ! reset the stabilizer for this depth
             self%Gsize(d+1) = 0
             self%loc(d+1) = -1
             self%loc(d) = self%loc(d) + 1
          end do
       end if

       ! If we are on the base layer of the top layer of the tree then
       ! we have created a lookup table that lets us know if a
       ! location on this layer of the tree is equivalent to a
       ! location we have already visited. If the next location we
       ! would visit is already equivalent to one we have already seen
       ! then we want to skip it. We only do this for the first layer
       ! in the tree because it allows us to eliminate the largest
       ! possible chunks of the tree possible and because implementing
       ! the method on all levels would result in a complex system of
       ! linked lookup tables between the levels of the tree. base is
       ! the lookup table for the first level in the tree and is
       ! initialized to be 0 if the base for the current location is
       ! anything else then we have seen it before and it needs to be
       ! skipped.
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

       ! Allocate the stabilizers for the next level to be the appropriate sizes.
       if (d > 0) then
          allocate(self%G%layer(d+1)%perms(self%Gsize(d),self%n),STAT=status)
          if(status/=0) stop "Allocation failed in increment_location: self%G%layer%perms."
          allocate(self%A%layer(d+1)%perms(self%Gsize(d),size(self%A%layer(d)%perms,2)),STAT=status)
          if(status/=0) stop "Allocation failed in increment_location: self%A%layer%perms."
          self%G%layer(d+1)%perms = 0
          self%A%layer(d+1)%perms = 0
       end if
    else
       ! There is only one atomic species in the system that we must
       ! consider so we have no where to go in the tree and it is
       ! fully explored.
       self%done = .True.
    endif
  endsubroutine increment_location

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
    !!<local name="perm_i">Indexing variable</local>
    !!<local name="perm_j">Indexing variable</local>
    !!<local name="stab_size">The size of the next stabilizer group.</local>
    !!<local name="new_loc">The location of the rotated branch.</local>
    !!<local name="min_stab">Temporary array for the stabilizers.</local>
    !!<local name="min_stab_a">Temporary array for arrow stabilizers.</local>
    !!<local name="d">Keeps track of how deep in the tree we are.</local>
    integer :: rotatedlabel(size(label))
    integer :: perm_i, perm_j, stab_size, new_loc, d, status
    integer, allocatable :: min_stab(:,:), min_stab_a(:,:)

    d = self%depth()
    ! Make sure that the number of stabilizers is initially 0
    self%Gsize(d + 1) = 0

    ! Now apply the symmetry group in order to determine if the
    ! labeling is unique and to find the stabilizers
    groupCheck: do perm_i = 1, self%Gsize(d)
       ! Apply the ith permutation to the labeling
       rotatedlabel = label(self%G%layer(d)%perms(perm_i,:))
       if (all( rotatedlabel == label )) then
          ! If we are at the labeling is the same then it is possible
          ! that this structure is superperiodic so we have to
          ! check. This will only be the case if we are on the bottom
          ! layer of the tree (d == k-1). We only want to check if
          ! this is a superperiodic structure if the perumtation is a
          ! translation of the lattice (i.e. it's one of the first n
          ! permutations where n is the number of unit cells in the
          ! system) and if the pemutation is not the identity
          ! (permutiation 1).
          if ((d == self%k -1) .and. (perm_i <= symsize) .and. (perm_i > 1) .and. (.not. fixedcell) .and. (self%nArrows == 0)) then
          ! if ((d == self%k -1) .and. (perm_i <= symsize) .and. (perm_i > 1) .and. (.not. fixedcell)) then
             ! If this permutation is the same after being acted on by
             ! one ofiu the first n permutations of the entire system
             ! then the labeling is superperiodic and should not be
             ! included in the list.
             do perm_j = 2, symsize
                if (all(self%G%layer(1)%perms(perm_j,:) == self%G%layer(d)%perms(perm_i,:))) then
                   self%unique = .False.
                   exit groupCheck
                end if
             end do
          end if
          ! Otherwise this permutation is a stabilizer, i.e., it fixes the
          ! current labeling
          stab_size = self%Gsize(d+1) + 1
          self%Gsize(d+1) = stab_size
          self%G%layer(d+1)%perms(stab_size,:) = self%G%layer(d)%perms(perm_i,:)
          self%A%layer(d+1)%perms(stab_size,:) = self%A%layer(d)%perms(perm_i,:)
       else 
          ! check to see if the configuration is unique by finding the
          ! location of the permuted labeling. If it comes before the
          ! original labeling then the original labeling is a
          ! duplicate of something we've already seen.
          call self%get_loc(rotatedlabel,new_loc)
          if (d == 1) then
             self%base(self%loc(d)+1) = 1
             if (new_loc < self%loc(d)) then
                ! This labeling appeared earlier in the list
                self%unique = .False.
             elseif (new_loc > self%loc(d)) then
                self%base(new_loc + 1) = 1
             end if
          else
             if (new_loc < self%loc(d)) then
                ! The rotated labeling is equivalent
                ! to an earlier configuration.
                self%unique = .False.
                exit groupCheck
             end if
          end if
       end if
    end do groupCheck
   
  END SUBROUTINE check_labeling

  !!<summary>This routine takes a labeling and returns it's location
  !!in the branch of the tree. This subroutine is an implementation of
  !!the algorithm described in section 3.1 of
  !!http://msg.byu.edu/papers/enum3.pdf (graphically displayed in Fig
  !!5).</summary>
  !!<parameter name="labeling" regular="true">The labeling whose
  !!location is to be determined.</parameter>
  !!<parameter name="index" regular="true">The output location.</parameter>
  SUBROUTINE generateLocationFromColoring(self,labeling,index)
    class(tree)          :: self
    integer, intent(in)  :: labeling(:)
    integer, intent(out) :: index

    !!<local name="new_labeling">Contains a binary labeling, i.e.,
    !!'1's and '0's.</local>
    !!<local name="site_i">Counting term index.</local>
    !!<local name="site_j">Counting term index.</local>
    !!<local name="nzeros">The number of zeros left in the labeling.</local>
    !!<local name="index">The output location.</local>
    !!<local name="ndigits">The number of digits to the current zero.</local>
    !!<local name="numones">The number of '1's to the right of the current zero.</local>
    !!<local name="minl">The location of the next zero.</local>
    !!<local name="lastone">The index of the last '1' in the binary array.</local>
    !!<local name="d">Keeps track of how deep in the tree we are.</local>
    integer, allocatable :: new_labeling(:)
    integer :: site_i, nzeros, numones, lastone, nLeft, d, status
    integer :: ndigits(1), minl(1)

    d = self%depth()
    ! The labeling needs to have only '1's for the current color, have
    ! all smaller integer colorings removed, and be zeros
    ! otherwise. We also keey track of the location of the last '1' in the list.
    nLeft = count(labeling == 0, 1) + count(labeling >= d, 1)
    allocate(new_labeling(nLeft),STAT=status)
    if(status/=0) stop "Allocation failed in generateLocationFromColoring: new_labeling."
    new_labeling = 0
    where(pack(labeling,labeling >= d .or. labeling == 0) == d)
       new_labeling = 1
    end where

    ! Use the binomials to find the index
    index = 0
    do site_i = 1, nLeft
       if (new_labeling(site_i) == 0) then
          index = index + binomial(nLeft-site_i, count(new_labeling(site_i:)==1)-1)
       end if
    end do
    deallocate(new_labeling)
  END SUBROUTINE  generateLocationFromColoring

  !!<summary>This subroutine takes a unique coloring and finds the
  !!unique arrow arrangements that exist for that coloring.</summary>
  !!<parameter name="coloring" regular="true">The integer array containing the
  !! labeling for the colors</parameter>
  !!<parameter name="symsize" regular="true">The size of the supercell.</parameter>
  !!<parameter name="nfound" regular="true">Counter for the total number
  !!of labelings.</parameter>
  !!<parameter name="scount" regular="true">Counter for the number of
  !!labelings for this size.</parameter>
  !!<parameter name="HNFcnt" regular="true">Counter for the HNFs.</parameter>
  !!<parameter name="iBlock" regular="true">Index in the permIndx
  !!corresponding to the current block of the HNFs.</parameter>
  !!<parameter name="hnf_degen" regular="true">The degeneracy of the HNFs.</parameter>
  !!<parameter name="fixOp" regular="true">Lattice fixing operations (type opList).</parameter>
  !!<parameter name="SNF" regular="true">List of the SNFs.</parameter>
  !!<parameter name="HNF" regular="true">List of the HNFs.</parameter>
  !!<parameter name="LT" regular="true">List of the left transforms.</parameter>
  !!<parameter name="permIndx" regular="true">List of the different permutations
  !!groups.</parameter>
  !!<parameter name="equivalencies" regular="true">The list of
  !!equivalencies of the system.</parameter>
  subroutine addArrowsToEnumeration(self,coloring,symsize,nfound,scount,HNFcnt,&
                        iBlock,hnf_degen,fixOp,SNF,HNF,LT,equivalencies,permIndx)
    class(tree) :: self
    integer, intent(in) :: coloring(:)
    integer, intent(in)      :: symsize, HNFcnt, iBlock
    integer, intent(inout)   :: scount, nfound
    integer, intent(in)      :: SNF(:,:,:), HNF(:,:,:), LT(:,:,:)
    integer, intent(in)      :: equivalencies(:), hnf_degen(:), permIndx(:)
    type(opList), intent(in) :: fixOp(:)    

    !!<local name="all_sites">The arrow label being tested, including
    !!places without arrows.</local>
    !!<local name="arrowing">The arrow label for locations with arrows</local>
    !!<local name="max_arrowings">An array of the maximum of each arrow
    !!allowed.</local>
    !!<local name="maxindex">The maximum index possible for the arrowings.</local>
    !!<local name="d">The in the tree.</local>
    !!<local name="tempcoloring">Temporary copy of the coloring.</local>
    !!<local name="rotated_arrowing">The full labeling after the site permutation
    !!has been applied</local>
    !!<local name="permuted_arrowing">The full labeling after the site permutation
    !!and arrow permutation have been applied</local>
    !!<local name="full_coloring">The labeling for the colors after the arrow species have
    !!been mapped back to the original ones.</local>
    !!<local name="newindex">The transmogrified arrowings index.</local>
    !!<local name="arrow_dim">The number of directions the arrows can point.</local>
    !!<local name="nperms">The total number of permutations for this level of the
    !!tree.</local>
    integer :: all_sites(size(coloring)), rotated_arrowing(size(coloring)), full_coloring(size(coloring))
    integer :: permuted_arrowing(size(coloring)), permuted_coloring(size(coloring))
    integer, allocatable :: max_arrowings(:), arrowing(:)
    integer :: maxindex, originalIndex, d, newindex
    integer :: perm_i, perm_j, site_i, site_j, status, arrow_dim, nperms

    d = self%depth()

    ! Find the number of directions by finding the largest number in
    ! the arrow perm.
    d = d + 1
    nperms = self%Gsize(d)
    arrow_dim = maxval(self%A%layer(d)%perms(1,:), 1)
    
    allocate(max_arrowings(self%nArrows),arrowing(self%nArrows),STAT=status)
    if(status/=0) stop "Allocation failed in addArrowsToEnumeration: max_arrowings, arrowings."
    ! Within this algorithm the arrowings are treated as an odometer
    ! that increments through all possible arrow arrangements for the
    ! system. Each arrowed site can have 0 to 5 for bulk enumerations
    ! or 3 for surface enumerations. The largest possible odemeter is
    ! then max_arrowings.
    max_arrowings = arrow_dim !-1
    call generateIndexFromArrowing(max_arrowings, arrow_dim, maxindex)
    ! Loop over every possible arrow array by going from the index of
    ! 0 to the maxindex found for the max_arrowings array.
    do originalIndex = 0, maxindex
       self%unique = .True.
       ! In order to apply the symmetry group we need to construct a
       ! labeling for every site even this that don't have arrows on them.
       all_sites = 0

       ! Get the arrowing that goes with this index
       call generateArrowingFromIndex(originalIndex, self%nArrows, arrow_dim, arrowing)
       ! Create a labeling that includes the sites that don't have arrows.
       site_j = 1
       do site_i = 1, size(coloring)
          if (any(self%color_map(:,1) == coloring(site_i)) .and. (self%nArrows >= site_j)) then
             all_sites(site_i) = arrowing(site_j) 
             site_j = site_j + 1
          end if
       end do

       ! Use the stabilizers to determine if the arrowing is unique.
       groupCheck: do perm_i = 1, nperms
          rotated_arrowing = all_sites(self%G%layer(d)%perms(perm_i,:))
          permuted_coloring = coloring(self%G%layer(d)%perms(perm_i,:))
          do site_i =1,size(permuted_arrowing)
             if (rotated_arrowing(site_i) /= 0) then
                permuted_arrowing(site_i) = self%A%layer(d)%perms(perm_i,rotated_arrowing(site_i))
             else
                permuted_arrowing(site_i) = 0                
             end if
          end do
          if (all(permuted_coloring == coloring) .and. all(permuted_arrowing == all_sites)) then
             if ((perm_i <= symsize) .and. (perm_i > 1)) then
                do perm_j = 2, symsize
                   if (all(self%G%layer(1)%perms(perm_j,:) == self%G%layer(d)%perms(perm_i,:)) &
                        .and. all(self%A%layer(1)%perms(perm_j,:) == &
                        self%A%layer(d)%perms(perm_i,:))) then
                      ! This arrangement is superperiodic and not unique.
                      self%unique = .False.
                      exit groupCheck
                   end if
                end do
             end if
          else
             call generateIndexFromArrowing(permuted_arrowing, arrow_dim, newindex)
             ! If the new index came before the original then this
             ! arrowing is a duplicate.

             if (newindex < originalIndex) then
                self%unique = .False.
                exit groupCheck
             end if
          end if
       end do groupCheck
       if (self%unique .eqv. .True.) then
          ! If the arrowing is unique then we need to undo our
          ! colormap so we can save the labeling to file.
          do site_i=1, size(coloring)
             if (any(self%color_map(:,1) == coloring(site_i))) then
                do site_j = 1, size(self%color_map,1)
                   if (self%color_map(site_j,1) == coloring(site_i)) then
                      full_coloring(site_i) = self%color_map(site_j,2) -1
                   end if
                end do
             else
                full_coloring(site_i) = coloring(site_i) -1
             end if
          end do
          call write_single_labeling(full_coloring,symsize,nfound,scount,HNFcnt,&
               iBlock,hnf_degen,fixOp,SNF,HNF,LT,equivalencies,permIndx,arrow_label=all_sites)

       end if
    end do
    
  end subroutine addArrowsToEnumeration

  !!<summary>Writes a single labeling to file for the enum4 code. This
  !!code is a duplicae of the write_labelings code except
  !!with the outermost do loop removed.</summary>
  !!<parameter name="n" regular="true">The size of the supercell.</parameter>
  !!<parameter name="Tcnt" regular="true">Counter for the total number
  !!of labelings.</parameter>
  !!<parameter name="Scnt" regular="true">Counter for the number of
  !!labelings for this size.</parameter>
  !!<parameter name="Hcnt" regular="true">Counter for the HNFs.</parameter>
  !!<parameter name="HNFi" regular="true">Index in the permIndx
  !!corresponding to the current block of the HNFs.</parameter>
  !!<parameter name="hnf_degen" regular="true">The degeneracy of the HNFs.</parameter>
  !!<parameter name="fixOp" regular="true">Lattice fixing operations (type opList).</parameter>
  !!<parameter name="SNFlist" regular="true">List of the SNFs.</parameter>
  !!<parameter name="HNFlist" regular="true">List of the HNFs.</parameter>
  !!<parameter name="L" regular="true">List of the left transforms.</parameter>
  !!<parameter name="permIndx" regular="true">List of the different permutations
  !!groups.</parameter>
  !!<parameter name="equivalencies" regular="true">The list of
  !!equivalencies of the system.</parameter>
  !!<parameter name="labeling" regular="true">The labeling to be written to file.</parameter>
  !!<parameter name="arrow_label" regular="true">(Optional) The arrow
  !!label to be written to file if present in the
  !!enumeration.</parameter>
  subroutine write_single_labeling(labeling,n,Tcnt,Scnt,Hcnt,HNFi,hnf_degen,fixOp,SNFlist,HNFlist,L,equivalencies,permIndx, arrow_label)
    integer, intent(in)      :: n, Hcnt, HNFi
    integer, intent(inout)   :: Scnt, Tcnt
    integer, intent(in)      :: SNFlist(:,:,:), HNFlist(:,:,:), L(:,:,:)
    integer, intent(in)      :: equivalencies(:), hnf_degen(:), permIndx(:), labeling(:)
    type(opList), intent(in) :: fixOp(:)
    integer, optional, intent(in):: arrow_label(:)
    
    !!<local name="conc_check">Are the concentrations restricted</local>
    !!<local name="lab_degen">The degeneracy of this label.</local>
    !!<local name="iHNF">Counter that loops over the HNFs</local>
    !!<local name="vsH">Vector for matching the HNF to the list.</local>
    !!<local name="struct_enum_out_formatstring">Format for output file.</local>
    !!<local name="dummy">Dummy string the keeps track if the size of the labeling.</local>
    !!<local name="status">Allocation exit flag</local>
    !!<local name="i">Loop variable</local>
    !!<local name="jHNF">Which HNF this is in the permIndx</local>
    !!<local name="nHNF">How many of this HNF are there.</local>
    logical :: conc_check = .true.
    integer :: lab_degen, iHNF, status, i, jHNF, nHNF
    integer, allocatable :: vsH(:)
    character(105) :: struct_enum_out_formatstring
    character(3) :: dummy

    ! Labeling Postprocessing data
    !!<local name="nALLD">Size of full dset.</local>
    !!<local name="allD">help array: (/1,2,3,....,nALLD/)</local>
    !!<local name="pplabeling">Labeling after postprocessing.</local>
    !!<local name="postprocessLabeling">True if we need to perform
    !!postprocessing.</local>
    !!<local name="allD2LabelD">Respective d-vector ID in the labeling dset.</local>
    integer              :: nAllD
    integer, allocatable :: allD(:)
    integer, allocatable :: pplabeling(:)
    logical :: postprocessLabeling
    integer, allocatable :: allD2LabelD(:)


    ! Postprocessing labelings: setup
    nAllD = size(equivalencies)
    allocate(allD(nAllD)); allD = (/ (i,i=1,nAllD) /)
    allocate(allD2LabelD(nAllD));
    ! 1) Check whether we have to postprocess the labeling before writing
    !    it out Postprocessing is needed if we do not want to enumerate
    !    all dset members of a primitive unit cell due to some
    !    equivalencies.  For example, in a 1x1 symmetric surface slab (
    !    (*) denotes an atom ):
    !
    !      ------------------------------------------------ surface
    !         (*)  topmost surface layer, dvector# 1
    !         (*)  dvector# 2
    !         (*)  dvector# 3
    !         (*)  bottommost surface layer, dvector# 4
    !      ------------------------------------------------ surface
    !
    !    the topmost and the bottommost atom should always have the same
    !    occupancy, as well as dvector 2 and dvector 3 should. You can
    !    therefore specify the following equivalency list:
    !
    !         equivalency of dvector# | 1 2 3 4
    !         ---------------------------------   
    !         equivalency list        | 1 2 2 1
    !
    !    which means that dvector# 1 and dvector# 4 have to have the same
    !    occupancy, they are equivalent by enumeration. The same is true
    !    for dvector# 2 and dvector# 3
    !
    !    In this example, the enumeration code should only find
    !    enumerations of dvector# 1 and dvector# 2.  The occupations of
    !    dvector# 3 and dvector# 4 are then constructed in a
    !    postprocessing step.  The postprocessing step takes the
    !    enumerated form (i.e. dvectors# 1 and 2, e.g. a labeling 0101 for
    !    two unit cells) and tranforms it into a form that is valid for
    !    ALL dvectors, e.g. labeling 01010101 for two unit cells).
    !
    postprocessLabeling = .not. (all(  abs( equivalencies-allD ) ==0))
    allocate(pplabeling(n*nAllD))
    ! this is ok whether we do postprocessing (nAllD>nD) or not
    ! (nAllD==nD) 2) if yes: - prepare the postprocessed labeling
    ! (pplabeling) - generate a map: full dset -> enum dset that tells
    ! me for a dset member of the full dset what dset member in the
    ! enumeration dset it corresponds to
    if (postprocessLabeling) then
       allD2LabelD = (/(count(pack(allD,allD<=equivalencies(i))==equivalencies),i=1,nAllD)/)
       ! this construction is best explained by an example: suppose we
       ! have a dset 1,2,3,4,5 and the equivalent list is 1,4,4,4,5
       ! (so, dset member 1,4,5 will be enumerated, having positions
       ! 1,2,3 in labelings) the map should accomplish: 1 -> 1, 2 ->
       ! 2, 3 -> 2, 4 -> 2, 5 -> 3 first, the pack operation take only
       ! those in the dset <= equivalent point, then, the number of
       ! truly enumerated points is counted (truly enumerated points
       ! are points for which dset=equivalencies
    endif

    lab_degen = 0
    nHNF = count(permIndx==HNFi)
    allocate(vsH(nHNF),STAT=status); if (status/=0) stop "Allocation failed in write_single_labelings: vsH"
    
    ! Packing...
    vsH = pack((/(i,i=1,size(HNFlist,3))/), HNFi==permIndx); 

    write(dummy,'(I3)') n*nAllD    

    if (present(arrow_label)) then
       struct_enum_out_formatstring = '(i11,1x,i9,1x,i7,1x,i8,1x,i8,1x,i11,1x,i3,2x,i4,2x,3(i2,1x),2x,6(i2,1x),2x,9(i4,1x),2x,'//trim(dummy)//'i1,4x,'//trim(dummy)//'i1)'
    else
        struct_enum_out_formatstring = '(i11,1x,i9,1x,i7,1x,i8,1x,i8,1x,i11,1x,i3,2x,i4,2x,3(i2,1x),2x,6(i2,1x),2x,9(i4,1x),2x,'//trim(dummy)//'i1)'
     end if

    if (postprocessLabeling) then
       ! see the comments at the beginning of the current routine
       call postprocess_labeling(n,nAllD,labeling,pplabeling,allD2LabelD) 
    else
       pplabeling = labeling ! nothing changes
    endif
    
    do iHNF = 1, nHNF ! Write this labeling for each corresponding HNF
       jHNF = vsH(iHNF) ! Index of matching HNFs
       Tcnt = Tcnt + 1; Scnt = Scnt + 1

       if (present(arrow_label)) then
          write(14,struct_enum_out_formatstring) &
               Tcnt, Hcnt+iHNF,hnf_degen(jHNF),lab_degen,lab_degen*hnf_degen(jHNF),&
               Scnt,n,size(fixOp(jHNF)%rot,3),SNFlist(1,1,jHNF),SNFlist(2,2,jHNF),&
               SNFlist(3,3,jHNF),HNFlist(1,1,jHNF),HNFlist(2,1,jHNF),HNFlist(2,2,jHNF),&
               HNFlist(3,1,jHNF),HNFlist(3,2,jHNF),HNFlist(3,3,jHNF),transpose(L(:,:,jHNF)),&
               pplabeling, arrow_label
       else
          write(14,struct_enum_out_formatstring) &
               Tcnt, Hcnt+iHNF,hnf_degen(jHNF),lab_degen,lab_degen*hnf_degen(jHNF),&
               Scnt,n,size(fixOp(jHNF)%rot,3),SNFlist(1,1,jHNF),SNFlist(2,2,jHNF),&
               SNFlist(3,3,jHNF),HNFlist(1,1,jHNF),HNFlist(2,1,jHNF),HNFlist(2,2,jHNF),&
               HNFlist(3,1,jHNF),HNFlist(3,2,jHNF),HNFlist(3,3,jHNF),transpose(L(:,:,jHNF)),&
               pplabeling
       end if
    enddo ! loop over HNFs

  contains

    !!<summary>Purpose: take an enumeration labeling (for selected,
    !!non-equivalent (by enumeration) dvectors) and construct the full
    !!labeling for all dvectors. It makes sure that two dvectors in
    !!the same primitive unit cell of the parent lattice get the SAME
    !!labeling always. See also comments at the beginning of
    !!write_labelings.</summary>
    !!<parameter name="nUC" regular="true">number of unit
    !!cells.</parameter>
    !!<parameter name="nAllD" regular="true">number of d-vectors in
    !!the new labeling.</parameter>
    !!<parameter name="oldlabeling" regular="true">the old
    !!labeling</parameter>
    !!<parameter name="newlabeling" regular="true">the wanna-be new
    !!labeling.</parameter>
    !!<parameter name="newD2oldD" regular="true">{ new D# }: a map
    !!newD -> oldD</parameter>
    subroutine postprocess_labeling(nUC,nAllD,oldlabeling,newlabeling,newD2oldD)
      integer, intent(in) :: nUC, nAllD          
      integer, intent(in) :: oldlabeling(:)      
      integer, intent(out):: newlabeling(:)      
      integer, intent(in) :: newD2oldD(:)        

      !!<local name="newlab_pos">The position in the new label.</local>
      !!<local name="oldlab_pos">The position in the old label.</local>
      !!<local name="iD">Variable for looping.</local>
      !!<local name="iUC">Variable for looping.</local>
      integer :: newlab_pos, oldlab_pos
      integer :: iD,iUC
      
      do iD=1,nAllD
         do iUC=1,nUC
            newlab_pos = (iD-1)*nUC + iUC               ! position in the new labeling
            oldlab_pos = (newD2oldD(iD)-1)*nUC + iUC    ! corresponding
            ! position in the old labeling
            newlabeling(newlab_pos) = oldlabeling(oldlab_pos)
         enddo
      enddo

    end subroutine postprocess_labeling

  end subroutine write_single_labeling

END MODULE tree_class
