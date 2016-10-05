module itertools
  implicit none
  private
  public vararray, vararray2d, cproduct

  !!<summary>Provides a structure for variable length arrays that need to have their
  !!cartesian product taken.</summary>
  !!<usage>
  !!type(vararray) single
  !!allocate(single)
  !!single%init(array)
  !!</usage>
  type vararray
     !!<member name="items">The array of items to take cartesian product over.</member>
     !!<member name="length">The number of items in @CREF[this.items].</member>
     integer, pointer :: items(:)
     integer :: length
  contains
    procedure, public :: init => vararray_init
  end type vararray

  !!<summary>Provides a structure for variable length 2-d arrays that need to have their
  !!cartesian product taken.</summary>
  !!<usage>
  !!type(vararray) single
  !!allocate(single)
  !!single%init(array)
  !!</usage>
  type vararray2d
     !!<member name="items">The array of items to take cartesian product over.</member>
     !!<member name="length">The number of items in @CREF[this.items].</member>
     integer, pointer :: items(:,:)
     integer :: length
  contains
    procedure, public :: init => vararray2d_init
    procedure, public :: finalize => vararray2d_finalize
  end type vararray2d
contains
  !!<summary>Initializes the array items and length property.</summary>
  subroutine vararray_init(self, array, length, alloc)
    class(vararray) :: self
    integer, target, optional, intent(in) :: array(:)
    integer, optional, intent(in) :: length
    logical, optional, intent(in) :: alloc

    logical :: nalloc
    !We need to see if we are *copying* the array, or just referencing it.
    if (present(alloc)) then
       nalloc = alloc
    else
       nalloc = .false.
    end if

    if (present(array)) then
       if (nalloc) then
          allocate(self%items(size(array, 1)))
          self%items = array
       else
          self%items => array
       end if
       self%length = size(self%items, 1)
    else
       allocate(self%items(length))
       self%length = length
    end if
  end subroutine vararray_init

  !!<summary>Initializes the array items and length property.</summary>
  subroutine vararray2d_init(self, array, alloc)
    class(vararray2d) :: self
    integer, target, intent(in) :: array(:,:)
    logical, optional, intent(in) :: alloc

    logical :: nalloc
    !We need to see if we are *copying* the array, or just referencing it.
    if (present(alloc)) then
       nalloc = alloc
    else
       nalloc = .false.
    end if

    if (nalloc) then
       allocate(self%items(size(array, 1), size(array,2)))
       self%items = array
    else
       self%items => array
    end if
    self%length = size(self%items, 1)
  end subroutine vararray2d_init

  !!<summary>Resets the vararray2d to be empty so that an existing instance can be re-used.</summary>
  subroutine vararray2d_finalize(self)
    class(vararray2d) :: self
    self%items => null()
    self%length = 0
  end subroutine vararray2d_finalize

  !!<summary>Builds the cartesian product of the specified arrays.</summary>
  !!<parameter name="elements" regular="true">A 1-D ragged-array of vararray instances to take the 
  !!cartesian product over.</parameter>
  !!<parameter name="presult" regular="true">A 2-D array with each row being a cartesian product entry
  !!with one item contributed from each of the vararrays in elements.
  !!</parameter>
  subroutine cproduct(elements, presult)
    class(vararray), allocatable, intent(in) :: elements(:)
    integer, allocatable, intent(inout) :: presult(:,:)

    integer :: i, prod 

    prod = 1
    do i=1, size(elements,1)
       prod = prod*elements(i)%length 
    end do

    !Now, we figure out what size presult will have. It will be a Nxlen(elements) array where
    !N=\prod_i len(elements_i)

    allocate(presult(prod, size(elements,1)))
    call rproduct(elements, presult, 1, 1, prod)
  end subroutine cproduct

  !!<summary>Recursively takes the product of the existing subsets in presult with the
  !!new subset specified by depth.</summary>
  !!<parameter name="subsets">A subset whose value need to be cartesian-multiplied with the existing
  !!subsets in presult.</parameter>
  !!<parameter name="presult">The total cproduct result to alter.</parameter>
  !!<parameter name="depth">Recursive depth on the function stack.</parameter>
  !!<parameter name="start">The row number in presult to start altering values for.</parameter>
  !!<parameter name="right">The product of the array sizes from this array to the right.</parameter>
  recursive subroutine rproduct(subsets, presult, depth, start, right)
    class(vararray), intent(in) :: subsets(:)
    integer, allocatable, intent(inout) :: presult(:,:)
    integer, intent(in) :: depth, start, right
    
    !!<local name="cursor">Points to the row in the presult that is currently being updated
    !!with the values of the subset being multiplied.</local>
    !!<local name="kids">The product of array lengths for subsets to the right of this one
    !!*not* including this one.</local>
    integer :: i, j, cursor
    integer :: kids
    integer :: items(subsets(depth)%length)
    items = subsets(depth)%items

    !First we need to get the product of the list sizes to the right of this one.
    cursor = start
    kids = right/subsets(depth)%length
    do i=1, subsets(depth)%length
       do j=1, kids
          presult(cursor+j-1,depth) = items(i)
       end do
       !Now we just fill in the spots to the left of this column in the presult table.
       if (depth < size(subsets, 1)) then
          call rproduct(subsets, presult, depth+1, cursor, kids)
       end if
       cursor = cursor + kids
    end do
  end subroutine rproduct
end module itertools
