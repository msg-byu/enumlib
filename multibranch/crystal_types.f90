MODULE crystal_types
use num_types
implicit none
private
public derivCryst, opList, RotPermList

type RotPermList
   integer(si), pointer :: perm(:,:) ! First index is the permutation number
   integer nL ! Number of permutations in the list (perm is nL x nAtoms)
   real(dp), pointer :: v(:) ! Lattice vectors that return the rotated d-vectors
                             ! back into the first cell
endtype RotPermList

type opList  ! A list of rotation operations
   real(dp), pointer :: rot(:,:,:)
   real(dp), pointer :: shift(:,:)
endtype opList

type derivCryst 
   integer :: diag(3)   ! diagonal elements of the SNF
   real(dp):: pLat(3,3) ! Parent lattice
   integer, pointer :: HNF(:,:)! Hermite Normal Forms for each SNF
   integer, pointer :: L(:,:,:)  ! Left SNF transformation
   integer, pointer :: label(:,:)! list of configurations
endtype derivCryst

ENDMODULE
