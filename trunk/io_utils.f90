MODULE io_utils
use num_types
use numerical_utilities
use vector_matrix_utilities
implicit none
private
public read_input
CONTAINS

!***************************************************************************************************
subroutine read_input(title,LatDim,pLV,nD,d,k,Nmin,Nmax,eps,full,label,digit)
character(80) :: title, pLatTyp, fullpart
integer,intent(out):: Nmin, Nmax, k, LatDim, nD
real(dp),intent(out) :: pLV(3,3), eps
real(dp), pointer :: d(:,:)
integer, pointer :: label(:,:), digit(:)

logical full, err
integer iD, i
character(100) line

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
read(10,*)  nD
allocate(d(3,nD),label(k,nD),digit(nD))
label = -1
! This next part is a bit messy but it makes the input file easy to set up 
! (no need for formatted reads from the file)
do iD = 1, nD ! loop over all the d-vectors
   call co_ca(10,err)
   read(10,'(a100)') line
   line = adjustl(line) ! Remove preceding blanks
   do i = 1,3 ! Loop over x,y,z coordinates of d-vector
      read(line,*) d(i,iD) ! Get a coordinate of the d-vector 
      line = adjustl(line(index(line," "):)) ! Throw away the number we just read in
   enddo

! Now read in the labels for this d-vector
   ! Make sure that there is at least one comment marker in the line
   line(100:100) = "#"
   ! Throw away the comment and append a "/" at the end of the remaining string
   line = trim(line(1:index(line,"#")-1))//"/"
!   print *,"starting string",line
   do i = 1, k ! Loop over the number of (possible) labels, exit when there are no more /'s
      if (index(line,"/")==0) &
        stop "The labels for each d-vectors should be formated as #/#/#... where 0<=#<k"
      read(line,*) label(i,iD)
      ! Sanity check on the input for the label (perhaps not sufficient but catches some errors)
      if (label(i,iD) > k-1 .or. label(i,iD) < 0) then
         write(*,'("Incorrect number for a label, ''",i2,"'', on d-vector #",i2)') label(i,iD), iD
         stop
      endif
      line = adjustl(line(index(line,"/")+1:)) ! remove the label that we just read in
      if(line=="") exit ! No more labels so go to the next d-vector
   enddo
   digit(iD) = i ! Store the number of labels that were specified for each d-vector
enddo
if(all(digit<k)) then
   write(*,'("digit: ",1x,80i2)') digit
   write(*,'(/"Number of labels in input file is insufficient for a ",i1,"-nary case")')k
   stop
endif
!do iD = 1,nD
!   write(*,'(10i2)') label(:,iD)
!enddo
!write(*,'("digits: ",10i2)') digit
!stop

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
