MODULE io_utils
use num_types
use enumeration_types
use numerical_utilities
use vector_matrix_utilities
use utilities_module, only: ucase
implicit none
private
public read_input, write_lattice_symmetry_ops, write_rotperms_list
CONTAINS

!***************************************************************************************************
subroutine read_input(title,LatDim,pLV,nD,d,k,eq,Nmin,Nmax,eps,full,label,digit)
character(80) :: title, pLatTyp, fullpart
integer,intent(out):: Nmin, Nmax, k, LatDim, nD
real(dp),intent(out) :: pLV(3,3), eps
real(dp), pointer :: d(:,:)
integer, pointer :: label(:,:), digit(:)
integer, pointer :: eq(:)

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
allocate(d(3,nD),label(k,nD),digit(nD),eq(nD))
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
! Should also check that no labels were repeated.
enddo
! Check that each label appears at least once
do i = 0,k-1
   if(all(label/=i))then
      write(*,'("Not all of the labels were used. Label ",i1," was never used")') i
   stop; endif
enddo
! Check that no label appears twice for one member of the dset
do iD = 1, nD
   do i = 0,K-1
      if(count(label(:,iD)==i)>1)then
         write(*,'("Label # ",i1," appears more than once for d-set # ",i2)') i,iD
      stop; endif
   enddo
enddo

call co_ca(10,err)
read(10,*) line
call ucase(line)
if (line(1:1)=='E') then
  call co_ca(10,err)
  read(10,*) eq(:)
else
  eq = (/(i,i=1,nD)/)
  backspace(10)
endif
call co_ca(10,err)
read(10,*) Nmin, Nmax
call co_ca(10,err)
read(10,*) eps
call co_ca(10,err)
read(10,*) fullpart
call ucase(pLatTyp)
call ucase(fullpart)
if (pLatTyp=='SURF') then; LatDim = 2
   if (.not. equal((/pLV(2,1),pLV(3,1)/),(/0._dp,0._dp/),eps)) &
        stop 'For "surf" setting, first component of second and third &
               & must be zero'
else if(pLatTyp=='BULK') then; LatDim = 3
else; stop 'Specify "surf" or "bulk" in input file';endif

if (fullpart=='FULL') then; full = .true.
else if(fullpart=='PART') then; full = .false.
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

!***************************************************************************************************
! This routine writes the symmetry operations of the multilattice (defined in lat.in to a file. This
! isn't used anywhere but might be nice for debugging and other checks.
SUBROUTINE write_lattice_symmetry_ops(rot,shift,key)
real(dp), intent(in) :: rot(:,:,:), shift(:,:)
character(2), optional :: key

character(2) twoDkey
integer iOp, nOp,i
twoDkey = "3D"
if (present(key)) then; if(key=="2D") twoDkey = "2D"; endif

if(twoDkey=="2D") then; open(11,file="symops_enum_2D.out",status="unknown"); else
open(11,file="symops_enum.out",status="unknown");endif
nOp = size(rot,3)
write(11,'("Number of symmetry operations: ",i2)') nOp
do iOp = 1,nOp
   write(11,'("Op #:",i2)') iOp
   write(11,'(3(f10.6,1x))') (rot(:,i,iOp),i=1,3)
   write(11,'("shift: ",3(f8.4,1x))') shift(:,iOp)
enddo
END SUBROUTINE write_lattice_symmetry_ops

!***************************************************************************************************
! Write out the information contained in a rotations-permutations list
SUBROUTINE write_rotperms_list(rpList,listfile)
type(RotPermList), intent(in) :: rpList
character(80), intent(in) :: listfile
integer nP, iP
open(11,file=listfile)
nP = rpList%nL
if(size(rpList%perm,1)/=rpList%nL) stop "rp list not initialized correctly (write_rotperms_list in io_utils)"
write(11,'("Number of permutations: ",i2)') nP
do iP = 1, nP
   write(11,'("Perm #: ",i2,1x,"Perm: ",40i2)') iP, rpList%perm(iP,:)
enddo
close(11)
END SUBROUTINE write_rotperms_list
END MODULE io_utils
