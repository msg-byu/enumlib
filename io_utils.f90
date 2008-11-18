MODULE io_utils
use num_types
use numerical_utilities
use vector_matrix_utilities
implicit none
private
public read_input, read_nth_struct_from_enumlist
CONTAINS

!***************************************************************************************************
SUBROUTINE read_nth_struct_from_enumlist(fname,title,bulksurf,k,p,HNF,labeling,eps)
character(*), intent(in) :: fname ! File to seach (assumes same format as struct_enum.out)
character(*), intent(out):: title, bulksurf, labeling
real(dp), intent(out) :: p(3,3), eps
real(dp), pointer :: dvec(:,:)
integer, intent(out) :: k, HNF(3,3)

integer :: nD, iline, strN, sizeN, nAt, pgOps, diag(3), a, b, c, d, e, f, L(3,3), ioerr, i
real(dp):: Binv(3,3)

open(11,file=fname,status='old',iostat=ioerr)
if(ioerr/=0)then; write(*,'("Input file doesn''t exist:",a80)') trim(fname);endif
! Read in the title from the struct_enum.out file
read(11,'(a80)') title; title = trim(title)

! Read in surf/bulk mode marker
read(11,'(a1)') bulksurf

! Read in the parent lattice vectors
do i = 1,3; read(11,*) p(:,i); enddo
call matrix_inverse(p,Binv)

! Read in the number of d-vectors, then read each one in
read(11,*) nD
allocate(dvec(3,nD))
do i = 1,nD; read(11,*) dvec(:,i); enddo

! Read in the number of labels, i.e., binary, ternary, etc.
read(11,'(i2)') k
read(11,*)
read(11,*) eps

! Skip to the specified structure number
do iline = 1, strN+1
   read(11,*)! dummy; print *,dummy
enddo
! Read in the info for the given structure
read(11,*) strN, sizeN, nAt, pgOps, diag, a,b,c,d,e,f, L, labeling
close(11)
!print *, strN, sizeN, nAt, pgOps, diag, a,b,c,d,e,f,L,labeling

! Define the full HNF matrix
HNF = 0
HNF(1,1) = a; HNF(2,1) = b; HNF(2,2) = c
HNF(3,1) = d; HNF(3,2) = e; HNF(3,3) = f

ENDSUBROUTINE read_nth_line_from_enumlist

!***************************************************************************************************
subroutine read_input(title,LatDim,pLV,nD,d,k,Nmin,Nmax,eps,full)
character(80) :: title, pLatTyp, fullpart
integer,intent(out):: Nmin, Nmax, k, LatDim, nD
real(dp),intent(out) :: pLV(3,3), eps
real(dp), pointer :: d(:,:)
logical full, err
integer iD

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
read(10,*)  nD
allocate(d(3,nD))
do iD = 1, nD
   call co_ca(10,err)
   read(10,*) d(:,iD)
enddo
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
