MODULE enumeration_utilities
use num_types
use vector_matrix_utilities
use numerical_utilities
implicit none
private
public  map_enumStr_to_real_space, cartesian2direct

CONTAINS

!***************************************************************************************************
SUBROUTINE map_enumStr_to_real_space(k, n, HNF, labeling, pLV, pBas, eps, sLV, aBas, spin, gIndx, x, L, S)
integer, intent(in) :: k ! Number of atom types (number of colors in the enumeration)
integer, intent(in) :: n ! Number of atoms in the unit cell of given structure
integer, intent(in) :: HNF(3,3) ! Hermite normal form corresponding to the superlattice
character(80), intent(in) :: labeling ! List, 0..k-1, of the atomic type at each site
real(dp), intent(in) :: pLV(3,3), eps ! parent lattice vectors (Cartesian coordinates), epsilon
real(dp), intent(in) :: pBas(:,:) ! 3xn array of the interior points of the parent lattice
real(dp), intent(out) :: sLV(3,3) ! Superlattice vectors (Cartesian coordinates)
real(dp), pointer :: aBas(:,:) ! Atomic positions (cartesian coordinates)
integer, pointer :: spin(:) ! Occupation variable of the positions
integer, pointer :: gIndx(:) ! Label ordinal for each atom (position in the labeling)
real(dp), pointer :: x(:) ! Concentration of each component
integer, intent(in):: L(3,3) ! HNF->SNF left transform. Need this to map r->G (Eq. 3 in first paper)
integer, intent(in):: S(3) ! Diagonal entries of the SNF

integer a, b, c, d, e, f ! elements of the HNF matrix
integer digit ! one of the labelings in the labeling
integer ic, z1, z2, z3 ! Counter over number of interior points, steps in the g table ("lattice coords")
!!real(dp):: HNFinv(3,3)!, gTab2Latt2Cart(3,3) ! Inverse of the HNF, convert gtable coords to Cartesian 
integer iAt, i, iD, nD
!real(dp) :: map2G(3,3) ! Transform for mapping a real vector into the group (r->G)
real(dp) :: greal(3)   ! Floating point representation of the group element components (g-vector)
integer  :: g(3)       ! Integer version of greal
!integer  :: gIndx(n*size(pBas,2)) ! (ordinal) index of g-vector in the group
!!real(dp) :: Ainv(3,3) ! Inverse of the parent lattice vectors
logical err

allocate(gIndx(n*size(pBas,2)))
gIndx=-1

nD = size(pBas,2)
! Define the non-zero elements of the HNF matrix
a = HNF(1,1); b = HNF(2,1); c = HNF(2,2)
d = HNF(3,1); e = HNF(3,2); f = HNF(3,3)
! Compute the superlattice vectors 
sLV = matmul(pLV,HNF)

! Find the coordinates of the basis atoms
allocate(aBas(3,n*nD))
allocate(spin(n*nD),gIndx(n*nD))

! Map HNF points to real space
!!call matrix_inverse(real(HNF,dp),HNFinv)

! Let's get the fattest basis (Minkowski reduction)
call reduce_to_shortest_basis(sLV,sLV,eps)

! Find each atomic position from the g-space information
ic = 0  ! Keep track of the number of points mapped so far
do iD = 1, nD ! Loop over the number at sites/parent cell (the d set)
do z1 = 0, a-1 ! For the limits on the loops, see the "interior_points.pdf" write-up
   do z2 = (b*z1)/a, c+(b*z1)/a - 1
      do z3 = z1*(d-(e*b)/c)/a+(e*z2)/c, f+z1*(d-(e*b)/c)/a+(e*z2)/c - 1
         ic = ic + 1 !; if (ic > n) stop "Problem in basis atoms..."
         ! Atomic basis vector in Cartesian coordinates
!!          aBas(:,ic) = matmul(HNFinv,(/z1,z2,z3/))+pBas(:,iD) 
         aBas(:,ic)=matmul(pLV,(/z1,z2,z3/))+pBas(:,iD)
!!         write(*,'("at #: ",i3,4x,"position: ",3(f7.3,1x))') ic,aBas(:,ic) 
         greal = matmul(L,(/z1,z2,z3/)) ! Map position into the group
         g = nint(greal) ! Convert the g-vector from real to integer
         if(.not. equal(greal,g,eps)) stop "map2G didn't work in map_enumStr_to_real_space"
         g = modulo(g,S) ! Bring the g-vector back into the first tile
         gIndx(ic) = (iD-1)*S(1)*S(2)*S(3)+g(1)*S(2)*S(3)+g(2)*S(3)+g(3)+1
      enddo
   enddo
enddo
enddo
if (ic /= n*nD) stop "ERROR: map_enumStr_to_real_space: Didn't find the correct # of basis atoms"

! Now map each position into the group so that the proper label can be applied
!write(*,'(3("L: ",3(i1,1x),/))') (L(i,:),i=1,3)
!write(*,'("SNF values: ",/,3(i1,1x))') S
!write(*,'("G indicies: ",/,10(i1,1x))') gIndx

allocate(x(k))
x = 0.0
if (mod(k,2)==0) then
   do iAt = 1, n*nD
!      print *,iAt,gIndx(iAt),labeling(gIndx(iAt):gIndx(iAt))
      !print *,ichar(labeling(gIndx(iAt):gIndx(iAt)))
      i = ichar(labeling(gIndx(iAt):gIndx(iAt)))-48
      digit = i-k/2 ! convert 0..k-1 label to spin variable -k/2..k/2
      x(i+1) = x(i+1) + 1  ! Keep track of the concentration of each atom type
      if (digit<0) then
         spin(iAt) = digit
      else
         spin(iAt) = digit+1 ! skip 0 as a spin if k is even
      endif
   enddo
else
   do iAt = 1, n*nD
      i = ichar(labeling(gIndx(iAt):gIndx(iAt)))-48
      spin(iAt) = i-k/2
      x(i+1) = x(i+1) + 1 ! Keep track of the concentration of each atom type
   enddo   
endif
x = x/real(n*nD,dp)
ENDSUBROUTINE map_enumStr_to_real_space

!***************************************************************************************************
subroutine cartesian2direct(sLV,aBas, eps) 
real(dp), intent(in)    :: sLV(3,3) ! Superlattice vectors (Cartesian coordinates)
real(dp), intent(inout) :: aBas(:,:) ! Atomic positions (cartesian coordinates first, then direct)
real(dp), intent(in)    :: eps

real(dp) :: sLVinv(3,3)
integer iAt, nAt

nAt = size(aBas,2)
call matrix_inverse(sLV,sLVinv)

!!! Convert aBas to DIRECT COORDINATES
do iAt=1,nAt
  aBas(:,iAt) = matmul(sLVinv,aBas(:,iAt)) ! Put positions into "direct" coordinates
  ! This keeps the atomic coordinates inside the first unit cell---
  ! not necessary but aesthetically pleasing.
  do while(any(aBas(:,iAt) >= 1.0_dp - eps) .or. any(aBas(:,iAt) < 0.0_dp - eps)) 
    aBas(:,iAt) = merge(aBas(:,iAt), aBas(:,iAt) - 1.0_dp, aBas(:,iAt) <  1.0_dp - eps) 
    aBas(:,iAt) = merge(aBas(:,iAt), aBas(:,iAt) + 1.0_dp, aBas(:,iAt) >= 0.0_dp - eps)
  enddo
end do
!!! End of conversion to DIRECT COORDINATES
end subroutine cartesian2direct

END MODULE enumeration_utilities
