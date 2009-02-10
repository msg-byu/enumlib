MODULE enumeration_types
use num_types
implicit none
private
public derivCryst, opList, LabelRotationList, RotPermList

type RotPermList
   integer, pointer :: perm(:,:) ! First index is the permutation number
   integer, pointer :: RotIndx(:) ! Which rotations in the list fix the superlattice
   integer nL ! Number of permutations in the list (perm is nL x nAtoms)
   real(dp), pointer :: v(:,:,:) ! Lattice vectors that return the rotated d-vectors
                                 ! back into the first cell 3xnDxnOp
endtype RotPermList

type opList  ! A list of rotation operations
   real(dp), pointer :: rot(:,:,:)
   real(dp), pointer :: shift(:,:)
endtype opList

type LabelRotationList ! A list of indices corresponding to label rotations
   integer, pointer :: lr(:)
endtype LabelRotationList

type cryst ! a structure with a "real-space" definition
   real(dp) :: LV(3,3)
   real(dp), pointer :: bas(:,:) ! 3xn array of atom positions
   integer, pointer :: aTyp(:) ! array of n integers indicating the atom type of each atomic basis
endtype cryst

! Not used yet. Not clear that it is useful.
type derivCryst 
   integer :: diag(3)   ! diagonal elements of the SNF
   real(dp):: pLat(3,3) ! Parent lattice
   integer, pointer :: HNF(:,:,:)! Hermite Normal Forms for each SNF
   integer, pointer :: L(:,:,:)  ! Left SNF transformation
   integer, pointer :: labeling(:,:)! list of configurations
endtype derivCryst

type derivStruct
   integer :: diag(3)   ! diagonal elements of the SNF
   real(dp):: pLat(3,3) ! Parent lattice
   real(dp), pointer :: dVec(:,:) ! a 3xn list of the parent basis vectors
   integer :: nD ! Number of dVectors
   integer :: HNF(3,3)! Hermite Normal Forms for each SNF
   integer :: L(3,3)  ! Left SNF transformation
   integer n ! index, number of parent cells in the super 
   integer, pointer :: labeling(:)! list of labels
   integer, pointer :: conc(:) ! Number of atoms of each type in labeling
endtype derivStruct
ENDMODULE
