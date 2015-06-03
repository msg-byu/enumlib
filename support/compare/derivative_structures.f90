!!<summary> This module encapsulates the data types and routines needed to compare lists of enumerated
!!structures. The basics of the enumeration algorithm are discussed Hart & Forcade PRB 77 224115 (26
!!June 2008). </summary>
MODULE derivative_structures
use num_types
use io_utils
implicit none
private
public enum_list

!!<summary> Stores a list of permutations. The columns are the individual elements of a permuations,
!!and the rows are the different permutations</summary>
type perm_list
   integer, pointer, dimension(:,:) :: l
endtype perm_list

!!<summary>This data type contains an index for the SNF, HNF, and L matrix of the structure, as well
!!as an index to which list of permutations apply to this structure. It also contains a labeling
!!(list of atoms, one atom for each site in the structure) that is unique to the structure. </summary>
type structure
!!<member name="SNF"> List of diagonal entries of the SNF matrix.First index is list position,
!!</member>
!!<member name="HNF"> List of lower triangular entries of the HNF matrix. First index is list position,
!!second index (1-6) is for the HNF entries.</member>
   integer, dimension(3):: SNF 
   integer, dimension(6):: HNF
!   type(perm_list), pointer, dimension(:) :: perm
   integer, pointer, dimension(:) :: labeling
endtype structure

!!<summary> Contains the list of rotations and fractional shifts that constitute the space group
!!operations for a lattice</summary>
type spaceGroup
   real(dp), pointer, dimension(:,:,:) :: rot
   real(dp), pointer, dimension(:,:)   :: shift
contains
  procedure, public :: initialize => spaceGroup_initialize
endtype spaceGroup

!!<summary>This type contains everything necessary to define a parent lattice to enumerate over. It
!!also contains the symmetries of the lattice. </summary>
type lattice
   !!<member name="LV"> Lattice vectors that define the unit cell of the lattice</member>
   !!<member name="dvec"> "d-vectors" defining the position of lattice points inside the unit
   !!cell. Should always have at least one member.</member>
   !!<member name="dlab"> List of atom types (integers) that are allowed on each d-vector site</member>
   !!<member name="SpcGrp"> The rotational symmetries of the lattice+basis defined by the rest of
   !!the members of this lattice type.</member>
   real(dp), dimension(3,3)           :: LV
   real(dp), pointer, dimension(:,:)  :: dvec
   type(dlabel), pointer, dimension(:):: dLab
   type(spaceGroup)                   :: SpcGrp
contains
  procedure, public :: initialize => lattice_initialize
end type lattice


!!<summary>This structure stores an enumerated list of derivative superstructures, with all of the
!!auxiliary information such as SNFs, HNFs, parent lattice, etc </summary>
type enum_list
   !!<member name="parLat"> lattice vectors, atomic sites in the unit cell (d-vectors), and the list
   !!of allowed labels (atom types) on each site.  </member>
   !!<member name="k"> The number of different atom types used in the
   !!enumeration. The number of "colors."</member>
   !!<member name="SNF">
   !!  <summary>The diagonal elements of Smith Normal Form matrix. The SNF characterizes
   !!the translational part of the group of symmetries that apply to the enumerated structures. Many
   !!different structures may have the same SNF. This element is a list of all of the SNFs that are
   !!needed in the enumerated list.</summary>
   !!  <dimension type="column" index="1">An integer index for a 3x3 SNF matrix.</dimension>
   !!  <dimension type="column" index="2,3">The 3x3 SNF matrix rows and columns.</dimension>
   !!</member>
   !!<member name="HNF"> <summary>Similar to the SNF, the Hermite Normal Form, applies to some subset of the
   !!enumerated list. An HNF describes the lattice vectors of a structure.</summary>
   !!  <dimension type="column" index="1">An integer index for a 3x3 HNF matrix.</dimension>
   !!  <dimension type="column" index="2,3">The 3x3 HNF matrix rows and columns.</dimension>
   !!</member>
   !!<member name="L"> 
   !!  <summary> The left transform matrix that derives a superstructure from parent lattice.</summary>
   !!  <dimension type="column" index="1">An integer index for a 3x3 left transform matrix.</dimension>
   !!  <dimension type="column" index="2,3">The 3x3 left transform matrix, rows and columns.</dimension>
   !!</member>  
   type(lattice) :: parLat
   integer :: k
   type(structure), pointer, dimension(:) :: str
contains
  procedure, public :: read_in_list
endtype enum_list


!!<summary>This type stores a list of integers that represents the kinds of atoms (colors) that are
!!allowed on the corresponding site.</summary>
type dlabel
   !!<member name="d">
   integer, pointer, dimension(:) :: d
endtype dlabel

CONTAINS

!!<summary>Reads in a file in *struct_enum.out* format and sets up the enumerated list</summary>
!!<comments>This was copied and pasted from the "find_match..." routine in
!!enumeration_utilities. There is also a similar couple of routines in io_utils. This implementation
!!should supersede those and they should me removed eventually.</comments>
subroutine read_in_list(this,fname)
class(enum_list) :: this
character(len=:), allocatable, intent(in) :: fname ! Name of file to read in from

!!integer nMin, nMax ! Enumeration runs over all cells of size nMin to max size of nMax
!!integer LatDim ! 2D or 3D parent lattice?
!!real(dp)  eps ! Finite precision parameter
!!real(dp) :: parLV(3,3) ! Lattice vector of the parent lattice
!!
!!logical fullLab  !?? 
!!integer, pointer :: cRange(:,:) ! concentration ranges
!!integer, pointer :: label(:,:)
!!integer, pointer :: digit(:)
!!integer, pointer :: equivalencies(:) ! List of atoms that are equivalent

real(dp) :: p(3,3) ! Lattice vectors of the parent lattice
integer ioerr ! Stores code for file i/o error
character(800) title, dummy  ! Comment in line 1 of the enum file, dummy variable
character(1) :: bulksurf ! Lattice type (bulk 3D or surf 2D)
integer nD ! Number of sites in the basis (i.e., number of points in the multilattice)
real(dp), pointer :: dvec(:,:) => null() ! Atomic basis vectors in the unit cell
integer k ! Number of colors/label types (i.e., binary, ternary, etc.)
real(dp)  eps ! Finite precision parameter
integer i, j, n ! generic loop counters
type(structure), pointer :: struct_list

open(13,file="readcheck_read_in_struct_enum_list.out")
write(13,'(5(/),"This file echoes the input for the read_in_list routine and")')
write(13,'("prints some of the computed quantities as well. It is intended to")')
write(13,'("be useful for debugging and troubleshooting formatting issues")')
write(13,'("***************************************************")')

open(11,file=fname,status='old',iostat=ioerr)
if(ioerr/=0)then; write(*,'("Input file doesn''t exist:",a)') trim(fname);stop;endif
! Read in the title from the struct_enum.out file
read(11,'(a80)') title; title = trim(title)
write(13,'("Title: ",a80)') title

! Read in surf/bulk mode marker
read(11,'(a1)') bulksurf
write(13,'("bulk/surf: ",a1)') bulksurf
!GHcall ucase(bulksurf)
!GHif((LatDim==2 .and. bulksurf=="B") .or. (LatDim==3 .and. bulksurf=="S")) then 
!GH   print *,"The bulk/surf parameter in "//trim(fname)//" not consistent"
!GH   write(*,'("Bulk/surf parameter is: ",a1)') bulksurf
!GH   write(*,'("LatDim input in call is: ",i1)') LatDim
!GH   stop
!GHendif

! Read in the parent lattice vectors
do i = 1,3; read(11,*) p(:,i); enddo
write(13,'("Parent Lattice Vectors:",/,3(f8.3,1x))') (p(i,:),i=1,3)
!if (.not. equal(matmul(Ainv,pLV),nint(matmul(Ainv,pLV)),eps)) stop "Parent lattices don't match"

! Read in the number of d-vectors, then read each one in
read(11,*) nD
write(13,'("Number of d-vectors: ",i3)') nD
!if(nD/=size(dset,2)) stop "Numbers of d-vectors don't match"
allocate(dvec(3,nD))
do i = 1,nD; read(11,*) dvec(:,i); enddo
write(13,'("d-vectors:",/,40(3f8.4,1x,/))') (dvec(:,i),i=1,nD) 


! Read in the number of labels, i.e., binary, ternary, etc.
read(11,'(i2)') k
write(13,'("Number of labels (k): ",i2)') k 

read(11,*) ! Skip starting and ending cell sizes
read(11,*) eps
write(13,'("Finite precision parameter (epsilon): ",g17.10)') eps 

do 
   read(11,*) dummy
   if(dummy=="start") exit
enddo
write(13,'("Found ""start"" label. Now reading in structures line-by-line ")')  
nS = 0 ! Count number of structures in the list
do 
   read(11,'(a800)',iostat=ioerr) dummy
   if (ioerr/=0) then ! end of file
      exit
   end if
   nS = nS + 1
enddo
!Now that we know how many structures are in the file, allocate the storage and back up in the file
!so that we can read them into memory.
do j = 1, nS
   backspace(11)
end do
allocate(struct_list,nS)




close(13)
close(11)
end subroutine read_in_list

!!<summary>Gets the space group for the parent lattice of an enum_list</summary>
!!<parameter name="LV"></parameter>
!!<parameter name="dvec"></parameter>
!!<parameter name="dLab"></parameter>
subroutine spaceGroup_initialize(this, LV, dvec, dLab)
  class(spaceGroup) :: this
  real(dp), dimension(3,3)           :: LV
  real(dp), pointer, dimension(:,:)  :: dvec
  type(dlabel), pointer, dimension(:):: dLab  
  
end subroutine spaceGroup_initialize

!!<summary>Reads the specified file to initialize a lattice. </summary>
subroutine lattice_initialize(this, filename)
  class(lattice) :: this
  character(*) :: filename
  !Read in lattice configuration from struct_enum.in
  
  call this%SpcGrp%initialize(this%LV, this%dvec, this%dLab)
end subroutine lattice_initialize

END MODULE derivative_structures
