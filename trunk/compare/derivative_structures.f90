!!<summary> This module encapsulates the data types and routines needed to compare lists of enumerated
!!structures. The basics of the enumeration algorithm are discussed Hart & Forcade PRB 77 224115 (26
!!June 2008). </summary>
MODULE derivative_structures
use num_types
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
   integer :: SNFidx
   integer :: HNFidx
   integer :: Lidx
   integer :: permIdx
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
   integer, pointer, dimension(:,:,:):: SNF
   integer, pointer, dimension(:,:,:):: HNF
   integer, pointer, dimension(:,:,:):: L
   type(perm_list), pointer, dimension(:) :: perm
   type(structure), pointer, dimension(:) :: str
endtype enum_list


!!<summary>This type stores a list of integers that represents the kinds of atoms (colors) that are
!!allowed on the corresponding site.</summary>
type dlabel
   !!<member name="d">
   integer, pointer, dimension(:) :: d
endtype dlabel

CONTAINS


!!<summary></summary>
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
