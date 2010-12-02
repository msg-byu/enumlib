program test_driver
use num_types
use enumeration_utilities
implicit none
character(80) fname, title, sfname, dummy
real(dp), dimension(3,3) :: pLV, sLV
integer,  dimension(3,3) :: SNF, L
real(dp), pointer :: aBas(:,:), dset(:,:)
integer, pointer :: aTyp(:), pLabel(:,:), HNF(:,:,:)
integer :: LatDim, match
real(dp) :: eps
eps = 1e-4

LatDim = 3
call getarg(1,dummy)
read(dummy,'(a80)') fname
call getarg(2,dummy)
read(dummy,'(a80)') sfname
if(iargc()/=2) stop "Need two arguments: poscar and struct_enum.out-type file"
!sfname = "struct_enum.out"
call read_poscar(fname,title,sLV,aBas,aTyp)
call get_HNF_of_derivative_structure_old(sfname,sLV,aBas,aTyp,pLV,dset,HNF,SNF,L,eps)
call get_gspace_representation(pLV,dset,sLV,aBas,aTyp,HNF,LatDim,pLabel,eps)
call find_match_in_structenumout(sfname,pLV,dset,HNF,SNF,LatDim,pLabel,match,eps) 
if(match/=0) then
write(*,'("Structure #: ",i8," is a match to the structure in ",a80)') match, adjustl(fname)
else
write(*,'("No match to the structure in ",a80)') adjustl(fname)

endif
endprogram test_driver

