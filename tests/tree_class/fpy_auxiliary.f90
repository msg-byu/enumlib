!!<fortpy version="1.7.5" codeversion="1.7.0" />
!!<summary>Auto-generated auxiliary module exposing interfaces to save
!!instances of user-derived types.</summary>
module fpy_auxiliary
  use fortpy
  use iso_c_binding, only: c_ptr, c_associated, c_loc, c_null_ptr

  implicit none
  private
  
  public auxsave, auxread

  !!<summary>Provides a single call interface for saving user-derived data types
  !!for single values and arrays.</summary>
  interface auxsave
     module procedure 
  end interface auxsave

  !!<summary>Provides a single call interface for reading user-derived data types
  !!for single values and arrays.</summary>
  interface auxread
     module procedure 
  end interface auxread

  !!<summary>Stores the prefix path and address of a pointer-type user-derived type
  !!variable that is being saved.</summary>
  type, public :: fpy_address
     !!<member name="prefix">The file prefix to the instances top-level variable.</member>
     !!<member name="address">64-bit adress of the pointer in memory.</member>
     character(len=:), allocatable :: prefix
     type(c_ptr) :: address
  end type fpy_address
contains




  !!<summary>Adds the specified value to the end of the integer-valued list.</summary>
  !!<parameter name="list">The integer-valued list to append to.</parameter>
  !!<parameter name="appendage">The extra value to append to the list.</parameter>
  subroutine append_address(list, appendage)
    type(fpy_address), allocatable, intent(inout) :: list(:)
    type(fpy_address), intent(in) :: appendage
    type(fpy_address), allocatable :: templist(:)

    allocate(templist(size(list,1)+1))
    templist(1:size(list,1)) = list
    templist(size(list,1)+1) = appendage

    call move_alloc(templist, list)
  end subroutine append_address

  !!<summary>Returns the integer index of the specified address in the stack.</summary>
  !!<parameter name="stack">Array of 64-bit addresses to pointers.</parameter>
  !!<parameter name="avalue">Address to search for in the stack.</parameter>
  integer function address_loc(stack, avalue, qprefix)
    type(fpy_address), intent(in) :: stack(:), avalue
    logical, optional, intent(in) :: qprefix
    integer :: i
    address_loc = -1
    do i=1, size(stack, 1)
       if ((.not.present(qprefix)) .and. c_associated(stack(i)%address, avalue%address)) then
          address_loc = i
          exit
       else if (present(qprefix) .and. stack(i)%prefix .eq. avalue%prefix) then
          address_loc = i
          exit
       end if       
    end do
  end function address_loc

  !!<parameter name="drange">The dimensionality of the variable at the current prefix context.</parameter>
  !!<parameter name="prefix">The file context (ending in -) of the variable to check data range on.</parameter>
  !!<parameter name="members">The name of a non-array variable that terminates the variable chain.</parameter>
  subroutine fpy_get_drange(drange, prefix, members)
    integer, allocatable, intent(out) :: drange(:)
    character(len=*), intent(in) :: prefix, members(:)

    integer :: D, i, j, k
    logical :: exists, memexists
    character(len=:), allocatable :: catstr
    character(50) :: istr
    
    !First, check the dimensionality of the range by checking for a member at the first position.
    exists = .true.
    memexists = .false.
    i = 0
    D = 0
    catstr = ''

    do while (i .lt. 7)
       if (len(catstr) .gt. 0) then
          catstr = catstr//'.1'
       else
          catstr = '1'
       end if

       memexists = .false.
       j = 1
       do while (j .le. size(members) .and. .not. memexists)
          inquire(file=prefix//catstr//'-'//trim(members(j)), exist=memexists)
          j = j+1
       end do
       
       i = i + 1
       if (memexists) then
          D = i
       end if
    end do
    if (D .eq. 0) return
    allocate(drange(D))

    do k=1, D
       i = 1
       exists = .true.
       do while(exists)
          write (istr, *) i
          catstr = ''
          
          do j=1, k
             if (len(catstr) .gt. 0) catstr = catstr//'.'
             
             if (j .lt. k) then
                catstr = catstr//'1'
             else
                catstr = catstr//trim(adjustl(istr))
             end if
          end do
          do j=k, D-1
             catstr = catstr//'.1'
          end do

          memexists = .false.
          j = 1
          do while (j .le. size(members) .and. .not. memexists)
             inquire(file=prefix//catstr//'-'//trim(members(j)), exist=memexists)
             j = j+1
          end do

          exists = memexists
          if (exists) then
             i = i+1
          else
             i = i-1
          end if
       end do
       drange(k) = i
    end do    
  end subroutine fpy_get_drange
  
  !!<summary>Returns the integer index of the address in @CREF[param.stack] that
  !!has the same prefix.</summary>
  integer function fpy_address_index(stack, prefix)
    type(fpy_address), intent(in) :: stack(:)
    character(len=*), intent(in) :: prefix
    integer :: i, ind
    character(len=:), allocatable :: tprefix

    fpy_address_index = 0
    ind = index(prefix, '/', .true.)
    if (ind .gt. 0) then
       tprefix = prefix(ind+1:)
    else
       tprefix = prefix
    end if

    do i=1, size(stack)
       if (stack(i)%prefix .eq. tprefix) then
          fpy_address_index = i
          exit
       end if
    end do
  end function fpy_address_index
  
  !!<summary>Saves the list of prefixes in order for deserializing later.</summary>
  !!<parameter name="stack">The stack to save to file.</parameter>
  !!<parameter name="folder">The path to the variable's saving folder.</parameter>
  subroutine fpy_save_addresses(stack, folder)
    type(fpy_address), intent(in) :: stack(:)
    character(len=*), intent(in) :: folder
    !!<local name="savelist">The list of prefixes to save, in order.</local>
    character(100) :: savelist(size(stack, 1))
    integer :: i

    do i=1, size(stack, 1)
       savelist(i) = stack(i)%prefix
    end do
    call pysave(savelist, folder//'.fpy.address')
  end subroutine fpy_save_addresses

  !!<summary>Restores the stack of fpy_addresses so that a variable's folder can be
  !!deserialized back via aux_read.</summary>
  subroutine fpy_read_address(stack, folder)
    type(fpy_address), intent(out), allocatable :: stack(:)
    character(len=*), intent(in) :: folder
    character(100), allocatable :: savelist(:)
    integer :: i
    logical :: fpy_success
    
    call fpy_read(folder//'.fpy.address', '#', savelist, fpy_success)
    if (.not. fpy_success) then
       allocate(stack(0))
    else
       allocate(stack(size(savelist)))
       do i=1, size(savelist)
          stack(i)%prefix = savelist(i)
          stack(i)%address = c_null_ptr
       end do
    end if
  end subroutine fpy_read_address
  
  !!<summary>Returns the 64-bit integer address of the pointer at @CREF[param.cloc].</summary>
  !!<parameter name="ploc">The object returned by calling the @CREF[loc] interface on
  !!the variable.</parameter>
  !!<parameter name="prefix">The file location prefix for the variable being located.</parameter>
  type(fpy_address) function fpy_get_address(ploc, prefix)
    type(c_ptr), intent(in) :: ploc
    character(len=*), intent(in) :: prefix
    integer :: ind

    fpy_get_address%address = ploc
    ind = index(prefix, '/', .true.)
    if (ind .gt. 0) then
       fpy_get_address%prefix = prefix(ind+1:)
    else
       fpy_get_address%prefix = "_"
    end if       
  end function fpy_get_address    
end module fpy_auxiliary
