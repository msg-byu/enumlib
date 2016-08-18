!!<summary>Defines a trees class to be used for enumerating derivative structures at a fixed
!!concentration.</summary>
MODULE tree_class
  use combinatorics
  use group_theory
  use enumeration_types
  use arrow_related
  implicit none
  private
  public tree


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
     !!<member name="narrows">The integer number of arrows in the
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
     integer          :: narrows = 0
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
  !!<parameter name="makeG" regular="true">Flag: True-> Generate group, False ->
  !!Generators already make a group</parameter>
  subroutine initializeTree(self,colors,generators, arrow_group, color_map, makeG)
    class(tree)         :: self
    integer, intent(in) :: colors(:)
    integer, intent(in) :: color_map(:,:)
    integer, pointer, intent(in) :: arrow_group(:,:)
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
    self%nArrows = size(color_map,1)
    allocate(self%color_map(size(color_map,1),size(color_map,2)))
    self%color_map = color_map
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
       allocate(self%G%layer(self%k),self%A%layer(self%k)) !Allocate number of perm
       !lists. Don't need one for last color but makes the code
       !cleaner. Perhaps rethink this later
       allocate(self%G%layer(1)%perms(size(generators,1),size(generators,2)))
       allocate(self%A%layer(1)%perms(size(arrow_group,1),size(arrow_group,2)))
       self%G%layer(1)%perms = generators
       self%A%layer(1)%perms = arrow_group
       
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
       allocate(self%A%layer(d+1)%perms(self%Gsize(d),size(self%A%layer(d)%perms,2))) 
       self%G%layer(d+1)%perms = 0
       self%A%layer(d+1)%perms = 0
    else ! We are at the bottom of the tree (2nd lowest layer, but
       ! lowest always has just one coloring)
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
          allocate(self%A%layer(d+1)%perms(self%Gsize(d),size(self%A%layer(d)%perms,2)))
          self%G%layer(d+1)%perms = 0
          self%A%layer(d+1)%perms = 0
       end if
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
    !!<local name="i">Indexing variable</local>
    !!<local name="j">Indexing variable</local>
    !!<local name="stab_size">The size of the next stabilizer group.</local>
    !!<local name="new_loc">The location of the rotated branch.</local>
    !!<local name="min_stab">Temporary array for the stabilizers.</local>
    !!<local name="min_stab_a">Temporary array for arrow stabilizers.</local>
    !!<local name="d">Keeps track of how deep in the tree we are.</local>
    integer :: rotatedlabel(size(label))
    integer :: i, j, stab_size, new_loc, d
    integer, allocatable :: min_stab(:,:), min_stab_a(:,:)

    d = self%depth()
    ! Make sure that the number of stibalizers is initially 0
    self%Gsize(d + 1) = 0

    ! Now apply the symmetry group in order to determine if the
    ! labeling is unique and to find the stabilizers
    groupCheck: do i = 1, self%Gsize(d)
       ! Apply the ith permutation to the labeling
       rotatedlabel = label(self%G%layer(d)%perms(i,:))
       if (all( rotatedlabel == label )) then
          if ((d == self%k -1).and. (i <= symsize) .and. (i > 1) .and. (.not. fixedcell)) then
             do j = 2, symsize
                if (all(self%G%layer(1)%perms(j,:) == self%G%layer(d)%perms(i,:))) then
                   self%unique = .False.
                   exit groupCheck
                end if
             end do
          end if
          ! This permutation is a stabilizer, i.e., it fixes the
          ! current labeling
          stab_size = self%Gsize(d+1) + 1
          self%Gsize(d+1) = stab_size
          self%G%layer(d+1)%perms(stab_size,:) = self%G%layer(d)%perms(i,:)
          self%A%layer(d+1)%perms(stab_size,:) = self%A%layer(d)%perms(i,:)
       else 
          ! check to see if the configuration is unique
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

    if (self%unique .eqv. .True.) then
       allocate(min_stab(stab_size,self%n))
       allocate(min_stab_a(stab_size,size(self%A%layer(1)%perms,2)))
       min_stab = self%G%layer(d+1)%perms(1:stab_size,:)
       min_stab_a = self%A%layer(d+1)%perms(1:stab_size,:)
       deallocate(self%G%layer(d+1)%perms)
       deallocate(self%A%layer(d+1)%perms)
       allocate(self%G%layer(d+1)%perms(stab_size,self%n))
       allocate(self%A%layer(d+1)%perms(stab_size,size(self%A%layer(d)%perms,2)))
       self%G%layer(d+1)%perms = min_stab
       self%A%layer(d+1)%perms = min_stab_a
       deallocate(min_stab,min_stab_a)
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
    !!<local name="d">Keeps track of how deep in the tree we are.</local>
    integer, allocatable :: new_labeling(:)
    integer :: i, j, nzeros, numones, lastone, nLeft, d 
    integer :: ndigits(1), minl(1)

    d = self%depth()
    ! The labeling needs to have only '1's for the current color, have
    ! all smaller integer colorings removed, and be zeros
    ! otherwise. We also keey track of the location of the last '1' in the list.
    nLeft = count(labeling .eq. 0,1) + count(labeling >= d,1)
    allocate(new_labeling(nLeft))
    new_labeling = 0
    where(pack(labeling,labeling >= d .or. labeling == 0) == d)
       new_labeling = 1
    end where

    ! Use the binomials to find the index
    index = 0
    do i = 1, nLeft
       if (new_labeling(i) == 0) then
          index = index + binomial(nLeft-i, count(new_labeling(i:)==1)-1)
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

    !!<local name="full_arrowing">The arrow label being tested, including
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
    integer :: full_arrowing(size(coloring)), rotated_arrowing(size(coloring)), full_coloring(size(coloring))
    integer :: permuted_arrowing(size(coloring))
    integer, allocatable :: max_arrowings(:), arrowing(:)
    integer :: maxindex, originalIndex, d, newindex
    integer :: i, j

    d = self%depth()
    
    allocate(max_arrowings(self%narrows),arrowing(self%narrows))
    max_arrowings = 5

    call generateIndexFromArrowing(max_arrowings,6,maxindex)

    ! Loop over every possible arrow array.
    do originalIndex = 0, maxindex
       self%unique = .True.
       full_arrowing = 1

       ! Get the arrowing that goes with this index
       call generateArrowingFromIndex(originalIndex,self%narrows,5,arrowing)

       ! Create a labeling that includes the sites that don't have arrows.
       j = 1
       do i = 1, size(coloring)
          if (any(self%color_map(:,1) == coloring(i))) then
             full_arrowing(i) = arrowing(j) + 1
             j = j + 1
          end if
       end do

       ! Use the stabilizers to determine if the arrowing is unique.
       groupCheck: do i = 1, self%Gsize(d)
          rotated_arrowing = full_arrowing(self%G%layer(d)%perms(i,:))
          permuted_arrowing = rotated_arrowing(self%A%layer(d)%perms(i,:))
          
          call generateIndexFromArrowing(permuted_arrowing -1,6,newindex)

          if (newindex < originalIndex) then
             self%unique = .False.
             exit groupCheck
          end if
       end do groupCheck
       if (self%unique .eqv. .True.) then
          do i=1, size(coloring)
             if (any(self%color_map(:,1) == coloring(i))) then
                full_coloring(i) = self%color_map(coloring(i)-size(self%color_map,1),2) -1
             else
                full_coloring(i) = coloring(i) -1
             end if
          end do

          call write_arrow_labeling(full_coloring,full_arrowing-1,symsize,nfound,scount,HNFcnt,&
               iBlock,hnf_degen,fixOp,SNF,HNF,LT,equivalencies,permIndx)

       end if
    end do
    
  end subroutine addArrowsToEnumeration

END MODULE tree_class
