PROGRAM HNF_counter
use num_types 
use enumeration_types
use derivative_structure_generator
use io_utils
use vector_matrix_utilities
use kgrid_utilities
implicit none

integer nMin, nMax ! Numbers of various things
integer k ! Number of colors/label types (i.e., binary, ternary, etc.)
integer LatDim ! 2D or 3D parent lattice?
integer nD ! Number of sites in the basis (i.e., number of points in the multilattice)
real(dp), pointer :: d(:,:) => null()
character(1) :: latTyp
real(dp)  eps, radius
real(dp) :: parLV(3,3), minkBasis(3,3) ! Parent lattice vectors and  mink reduced basis
real(dp) :: kmesh(3,3)


character(80) title, fname
logical fullLab,concCheck
integer, pointer :: cRange(:,:)
integer, pointer :: label(:,:)
integer, pointer :: digit(:)
integer, pointer :: equivalencies(:)

integer ivol, i, iHNF
integer, pointer, dimension(:,:,:) :: HNF => null()
type(RotPermList) :: ParRPList ! Just the list for the parent lattice
integer, pointer, dimension(:,:,:) :: rdHNF =>null()
type(opList), pointer :: fixOp(:)  ! Symmetry operations that leave a multilattice unchanged
type(RotPermList), pointer :: RPList(:) ! Master list of the rotation permutation lists
!type(RotPermList), pointer :: rdRPList(:) ! Master list of the *unique* rotation permutation lists
real(dp), pointer :: uqlatts(:,:,:) => null()
integer, pointer :: hnf_degen(:)
!real(dp), pointer :: pd(:,:) ! d-vectors defining the multilattice (primary lattice only)


if (iargc()>=1) then
   call getarg(1,fname)
else
   fname = "struct_enum.in"
endif
call read_input(title,LatDim,parLV,nD,d,k,equivalencies,nMin,nMax,eps&
     &,fullLab,label,digit,fname,cRange,concCheck) ! Read in parent lattice vectors, etc.
if (LatDim==3) then; latTyp='b';else;latTyp='s';endif
!call gen_multilattice_derivatives(title, parLV,nD,d,k,nMin,nMax,latTyp,eps,fullLab,&
!         label,digit,equivalencies,concCheck,cRange)
call get_dvector_permutations(parLV,d,ParRPList,LatDim,eps)

open(14,file="HNF_list.out")
do ivol = nMin, nMax
!   write(*,'("Index: ",i4)') ivol
!   write(14,'("Index: ",i4)') ivol
   if (LatDim==3) then !<<< 2D? or 3D? >>
      call get_all_HNFs(ivol,HNF)    ! 3D
   else
      call get_all_2D_HNFs(ivol,HNF) ! 2D
   endif
   ! The *inverse* of each HNF, multiplied by the "parent" LVs (the
   ! reciprocal vectors) will define a kgrid.
   call remove_duplicate_kgrids(parLV,HNF,rdHNF,uqlatts,eps)

   !call remove_duplicate_lattices(HNF,LatDim,parLV,d,ParRPList,rdHNF,fixOp,RPList,uqlatts,hnf_degen,eps)
!   write(*,'("# of HNFs",i5)') size(rdHNF,3)
   ! Now we have a list of lattices (uqlatts). If we mink-reduce each
   ! one, then one of the new basis vectors will be the shortest basis
   ! vector for the lattice. 1/2 this distance will the radius of the
   ! largest touching spheres that we can fit in the lattice.
   do iHNF = 1, size(uqlatts,3)
      call minkowski_reduce_basis(uqlatts(:,:,iHNF),minkBasis,1e-10_dp)
!      print *,determinant(uqlatts(:,:,iHNF))
      radius = norm(minkBasis(:,1))
      if (norm(minkBasis(:,2)) < radius) radius = norm(minkBasis(:,2))
      if (norm(minkBasis(:,3)) < radius) radius = norm(minkBasis(:,3))
      radius = radius/2 ! Radius of the largest "touching" spheres
      ! that fit into the lattice
!      print *,radius
      write(14,'(i3,2x,i5,2x,f8.3,3x,6(1x,i3),3x,3(3(7.5,1x),1x))') ivol, iHNF,4*3.1415926/3*radius**3&
           &/abs(determinant(uqlatts(:,:,iHNF))),(/rdHNF(1,1&
           &,iHNF),rdHNF(2,1:2,iHNF),rdHNF(3,:,iHNF)/), (uqlatts(:,i&
           &,iHNF),i=1,3)
!      write(14,'(3(i3,1x))') (rdHNF(i,:,iHNF),i=1,3)
!      write(14,*)
!      write(14,'(3(f5.2,1x))') (uqlatts(i,:,iHNF),i=1,3)
!      write(14,*)
   end do
      
end do
close(14)

END PROGRAM HNF_counter
