MODULE io_utils
use num_types
use enumeration_types
use numerical_utilities
use vector_matrix_utilities
use utilities_module, only: ucase, ralloc
use rational_mathematics, only: HermiteNormalForm, SmithNormalForm
implicit none
private
public read_input, write_lattice_symmetry_ops, write_rotperms_list, read_in_cells_from_file
CONTAINS

subroutine read_in_cells_from_file(n,HNFList,pLat,eps)
integer, intent(in) :: n ! current index of superlattices
integer, pointer    :: HNFList(:,:,:) ! Output
real(dp), intent(in):: pLat(3,3), eps

integer status, ns, is, iv, i, cStr
logical err
character(80) line
real(dp), dimension(3,3)::  H, pLatInv
integer,  dimension(3,3)::  Hout, Hin, T, L, R, S
real(dp), allocatable :: inStrs(:,:,:)

ns = 0
open(43,file="fixed_cells.in")
open(44,file="debug_read_in_cells_from_file.out")
status=0 
! Count the number of structures in the input file
do while (status==0)
   call co_ca(43,err)
!   if(err) stop "Trouble reading in structures in read_in_cells_from_file"
   read(43,*,iostat=status) line
   if(index(line,"<str>")>0) ns = ns + 1
enddo
rewind(43)
call matrix_inverse(pLat,pLatInv,err)
if(err) stop "Matrix inversion failed in read_in_cell_from_file"
write(44,'("Inverse matrix of the parent lattice (vectors in columns):")')
do i = 1,3
   write(44,'(3(f12.6,1x))') pLatInv(i,:)
enddo

! Now read in the unit cells of the structures
allocate(HNFList(3,3,ns),inStrs(3,3,ns))
cStr = 0 ! Number of structures of the correct volume factor
do is = 1, ns
   call co_ca(43,err)
   i = 0
   do
      i = i + 1
      read(43,*,iostat=status) line
      if (status/=0 .or. i >100) stop " Didn't find an <str> flag in fixed_cells.in"
      if(index(line,"<str>")>0) exit
   enddo
   do iv = 1, 3
      call co_ca(43,err)
      read(43,*) inStrs(:,iv,is)
   enddo
   write(44,'(/,i2,"-th read-in structure (vectors in columns):")') is
   do i = 1,3
      write(44,'(3(f12.6,1x))') inStrs(i,:,is)
   enddo

   ! Now that we have the lattice vectors, let's check to see that
   ! the read-in structure is a derivative of the parent lattice
   H = matmul(pLatInv,inStrs(:,:,is))
   if(.not. equal(H,nint(H),eps)) then 
      write(*,'("The ",i2,"-th structure in the fixed_cells.in file is &
      & not a derivative structure")') is
      stop
   endif
   Hin = nint(H)
   write(44,'("Integer transform to generate the superlattice:")')
   do i = 1,3
      write(44,'(3(i2,1x))') Hin(i,:)
   enddo

   call HermiteNormalForm(Hin,Hout,T)
   write(44,'("HNF to generate the superlattice:")')
   do i = 1,3
      write(44,'(3(i2,1x))') Hout(i,:)
   enddo

   call SmithNormalForm(Hout,L,S,R)
   write(44,'("Smith Normal Form of the HNF:")')
   do i = 1,3
      write(44,'(3(i2,1x))') S(i,:)
   enddo
   write(44,'("Index of the superlattice:",i3)') S(1,1)*S(2,2)*S(3,3)
   if (S(1,1)*S(2,2)*S(3,3)==n) then
      cStr = cStr + 1
      HNFList(:,:,cStr) = Hout
   endif
enddo
HNFList => ralloc(HNFList,cStr)
close(43); close(44)
endsubroutine read_in_cells_from_file
!***************************************************************************************************
subroutine read_input(title,LatDim,pLV,nD,d,k,eq,Nmin,Nmax,eps,full&
     &,label,digit,fname,cRange,conc_check)
character(80) :: title, pLatTyp, fullpart
character(80), optional :: fname
integer,intent(out):: Nmin, Nmax, k, LatDim, nD
real(dp),intent(out) :: pLV(3,3), eps
real(dp), pointer :: d(:,:)
integer, pointer :: label(:,:), digit(:)
integer, pointer :: eq(:), cRange(:,:)
logical, intent(out):: conc_check

logical full, err
integer iD, i, status
character(100) line

open(99,file="readcheck_enum.out")
if(.not. present(fname)) then
   open(10,file='struct_enum.in',status='old')
else
   open(10,file=fname,status='old')
endif
call co_ca(10,err)
read(10,'(a80)') title
write(99,'("Title: ",a80)') title
call co_ca(10,err)
read(10,'(a4)') pLatTyp
call ucase(pLatTyp)
write(99,'("Lattice type (bulk or surface): ",a4)') pLatTyp
call co_ca(10,err)
read(10,*) pLV(:,1)
call co_ca(10,err)
read(10,*) pLV(:,2)
call co_ca(10,err)
read(10,*) pLV(:,3)
write(99,'(3(f7.3,1x))') transpose(pLV)
call co_ca(10,err)
read(10,*) k
write(99,'("Number of labels: ",i2)') k
call co_ca(10,err)
read(10,*)  nD
write(99,'("Number of d-vectors: ",i3)') nD
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
   write(99,'(3(f8.4,1x))',advance="no") d(:,iD)
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
      write(99,'(i1,"/")',advance="no") label(i,iD)
      line = adjustl(line(index(line,"/")+1:)) ! remove the label that we just read in
      if(line=="") exit ! No more labels so go to the next d-vector
   enddo
   write(99,*)
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
write(99,'("Min and Max cell sizes: ",2(i2,1x))') Nmin, Nmax
call co_ca(10,err)
read(10,*) eps
write(99,'("Epsilon: ",g12.4)') eps 
call co_ca(10,err)
read(10,*) fullpart
fullpart = adjustl(fullpart)
call ucase(fullpart)
write(99,'("full/part mode: ",a4)') fullpart

if (pLatTyp(1:4).eq.'SURF') then; LatDim = 2
   if (.not. equal((/pLV(2,1),pLV(3,1)/),(/0._dp,0._dp/),eps)) &
        stop 'For "surf" setting, first component of second and third &
               & must be zero'
else if(pLatTyp(1:4).eq.'BULK') then; LatDim = 3
else; stop 'Specify "surf" or "bulk" in input file';endif

if (fullpart(1:4).eq.'FULL') then; full = .true.
else if(fullpart(1:4).eq.'PART') then; full = .false.
else; stop 'Specify "full" or "part" in the input file';endif

! Read in the concentration ranges
allocate(cRange(k,3))
cRange = 0
do i = 1, k
   call co_ca(10,err)
   read(10,*,iostat=status) cRange(i,:)
   conc_check = .true.
   if (status/=0) then ! concentration is not specificed
      write(*,'("Concentration ranges are not specified")')
      cRange = 0 
      conc_check = .false.
      if (i>1) then
         write(*,'(//,"--<< WARNING: Concentration ranges are partially specified   >>--")')
         write(*,'(   "--<< WARNING: If you intended to specify them, fix the input >>--",//)')
      endif
      exit
   endif
enddo
open(98,file="debug_conc_check.out")
close(98,status="delete")
open(98,file="debug_site_restrictions.out")
close(98,status="delete")


! Write to the debug file
if (conc_check) then
   write(99,'("Concentration ranges are specified. Will run with using &
   & the fixed-concentration algorithm.")')
   do i = 1, k
      write(99,'("Type",i2," conc:",i2,"/",i2,"--",i2,"/",i2)') i,cRange(i,(/1,3/)),cRange(i,2:3)
   enddo
else
   write(99,'("Concentration ranges are *not* specified. Using &
   & the full-conc-range algorithm (original enum algorithm).")')
endif
close(10)

if (any(cRange<0)) stop "ERROR: negative input on concentrations in read_input"
do i = 1, k
   if (maxval(cRange(i,1:2))>cRange(i,3)) then
      write(*,'("ERROR: Numerator is larger than denominator.")')
      write(*,'("ERROR: Check the concentration input for element #:",i2)') i
      stop
   endif
enddo

write(99,'(/,"--<< Successfully read the struct_enum.in file >>--",/)')
close(99)
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

if(twoDkey=="2D") then; open(11,file="symops_enum_2D_parent_lattice.out",status="unknown"); else
open(11,file="symops_enum_parent_lattice.out",status="unknown");endif
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
write(11,'("Number of permutations: ",i4)') nP
do iP = 1, nP
   write(11,'("Perm #: ",i5,1x,"Perm: ",40(i4,1x))') iP, rpList%perm(iP,:)
enddo
close(11)
END SUBROUTINE write_rotperms_list
END MODULE io_utils
