MODULE enumeration_types
use num_types
implicit none
private
public derivCryst, opList, LabelRotationList

type opList  ! A list of rotation operations
   real(dp), pointer :: rot(:,:,:)
endtype opList

type LabelRotationList ! A list of indices corresponding to label rotations
   integer, pointer :: lr(:)
endtype LabelRotationList

type derivCryst 
   integer :: diag(3)   ! diagonal elements of the SNF
   real(dp):: pLat(3,3) ! Parent lattice
   integer, pointer :: HNF(:,:)! Hermite Normal Forms for each SNF
   integer, pointer :: L(:,:,:)  ! Left SNF transformation
   integer, pointer :: label(:,:)! list of configurations
endtype derivCryst

ENDMODULE
