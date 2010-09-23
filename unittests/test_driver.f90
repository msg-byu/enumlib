#if 0
!----------------------------------------------------------
! This is a test driver which can impliment several modules.
! 
! This will test the individual functions and give output where code is
! broken for debugging purposes.
!
! Author: Derek Carr
! Email: zaphinath@gmail.com
!
!----------------------------------------------------------
#END IF

PROGRAM test_driver
use test_driver_dsg
!use derivative_structure_gen_test


implicit none

integer      :: istat, mode
character(2) :: modein
character(1) :: overridein
character(100), dimension(:), allocatable :: cmode
integer :: imode, jmode, nmode
integer :: iRef

character(80) :: filename, header(4)

!--------- Test Menu ---------
nmode = 50; ! Number of modules
allocate (cmode(nmode))

cmode = ""
cmode(01) = "Test derivative_structure_generator.f90"
cmode(02) = "Test something new"

!--------- Write Menu ----------
write(*,'(A)')	"#"
write(*,'(A)')	"# UNIT TEST 4 UNCLE"
write(*,'(A)')	"#"
write(*,'(A)')	"# Menu:"

imode=0
do 
  imode=imode+1
  if (imode>nmode) then; write(*,'(A)') "#"; write(*,'(A)') ""; exit; endif ! no more modes? exit!
  if (cmode(imode)=="") then ! empty mode, i.e. no menu entry? print 1 blank line and go on with the next entry
    do jmode=imode,nmode
      if (cmode(jmode)/="") then; write(*,'(A)') "#"; imode=jmode-1; exit; endif
    enddo
  else
    write(*,'(2A)') "#   ", cmode(imode) ! print menu entry
  endif
enddo
!endif 

!TEST_override=.false.
if (iargc()>=1) then
  call getarg(1,modein)
  read(modein,*) mode
  !if (iargc()==2) then
    !call getarg(2,overridein)
    !if (overridein=="!") TEST_override=.true.
  !endif
endif


select case(mode)
case(01)
	write(*,'(A)')	"Case 01"
	call test_all_dsg()

case(02)
	write(*,'(A)')	"Case 02"
	
case default
	!case(01)
	!case(02)
end select

END PROGRAM test_driver


