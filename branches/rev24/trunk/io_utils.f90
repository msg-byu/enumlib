MODULE io_utils
use num_types
use numerical_utilities
use vector_matrix_utilities
implicit none
private
public read_input

CONTAINS
!***************************************************************************************************
subroutine read_input(title,LatDim,pLV,k,Nmin,Nmax,eps,full)
character(80) :: title, pLatTyp, fullpart
integer,intent(out):: Nmin, Nmax, k, LatDim
real(dp),intent(out) :: pLV(3,3), eps
logical full
logical err
open(10,file='struct_enum.in',status='old')
call co_ca(10,err)
read(10,'(a80)') title
call co_ca(10,err)
read(10,'(a4)') pLatTyp
call co_ca(10,err)
read(10,*) pLV(:,1)
call co_ca(10,err)
read(10,*) pLV(:,2)
call co_ca(10,err)
read(10,*) pLV(:,3)
call co_ca(10,err)
read(10,*) k
call co_ca(10,err)
read(10,*) Nmin, Nmax
call co_ca(10,err)
read(10,*) eps
call co_ca(10,err)
read(10,*) fullpart
if (pLatTyp=='surf') then; LatDim = 2
   if (.not. equal(pLV(:,1),(/1._dp,0._dp,0._dp/),eps)) &
        stop 'For "surf" setting, first vector must be 1,0,0'
   if (.not. equal((/pLV(2,1),pLV(3,1)/),(/0._dp,0._dp/),eps)) &
        stop 'For "surf" setting, first component of second and third &
               & must be zero'
else if(pLatTyp=='bulk') then; LatDim = 3
else; stop 'Specify "surf" or "bulk" in input file';endif

if (fullpart=='full') then; full = .true.
else if(fullpart=='part') then; full = .false.
else; stop 'Specify "full" or "part" on line 9 of input file';endif

end subroutine read_input

subroutine co_ca(unit,error)
!           subroutine was taken from the code of Ralf Drautz
!
!           co_ca cares about comments and blanks and avoids reading them,
!	        comment lines start with a #
implicit none
character(50) phrase !letter: contains the first letter of every line
logical   com !true if comment line is found
integer  unit, i, ios !unit specifies the unit to read from
logical error

com = .true.; error = .false.
do while ( com .eqv. .true.)
   read(unit,50,iostat=ios) phrase
   if (ios/=0) exit
   !              blank line ?
   if (phrase .ne. ' ') then
      i = index(phrase, '#')
      !                 # not found ?
      if (i .ne. 0) then
         !                    # first letter in line ?
         if (i .ne. 1) then
            phrase = phrase(1:i-1)
            if (phrase .ne. ' ') com= .false.
         endif
      else
         com = .false.
      endif
   endif
end do
50 format(50a)
backspace unit
if (ios/=0) error = .true.
end subroutine co_ca

END MODULE io_utils
