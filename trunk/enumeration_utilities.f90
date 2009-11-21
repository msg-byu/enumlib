MODULE enumeration_utilities
use num_types
use vector_matrix_utilities
use numerical_utilities
implicit none
private
public  map_enumStr_to_real_space
CONTAINS
!***************************************************************************************************
SUBROUTINE map_enumStr_to_real_space(k, n, HNF, labeling, pLV, pBas, eps, sLV, aBas, spin, x, L, S)
integer, intent(in) :: k ! Number of atom types (number of colors in the enumeration)
integer, intent(in) :: n ! Number of atoms in the unit cell of given structure
integer, intent(in) :: HNF(3,3) ! Hermite normal form corresponding to the superlattice
character(80), intent(in) :: labeling ! List, 0..k-1, of the atomic type at each site
real(dp), intent(in) :: pLV(3,3), eps ! parent lattice vectors (Cartesian coordinates), epsilon
real(dp), intent(in) :: pBas(:,:) ! 3xn array of the interior points of the parent lattice
real(dp), intent(out) :: sLV(3,3) ! Superlattice vectors (Cartesian coordinates)
real(dp), pointer :: aBas(:,:) ! Atomic positions (Cartesian coordinates)
integer, pointer :: spin(:) ! Occupation variable of the positions
real(dp), pointer :: x(:) ! Concentration of each component
integer, intent(in):: L(3,3) ! HNF->SNF left transform. Need this to map r->G (Eq. 3 in first paper)
integer, intent(in):: S(3) ! Diagonal entries of the SNF

integer a, b, c, d, e, f ! elements of the HNF matrix
integer digit ! one of the labelings in the labeling
integer ic, z1, z2, z3 ! Counter over number of interior points, steps in the g table ("lattice coords")
real(dp):: Sinv(3,3), gTab2Latt2Cart(3,3) ! Inverse of the superlattice, convert gtable coords to Cartesian 
integer iAt, i, iD, nD
real(dp) :: map2G(3,3) ! Transform for mapping a real vector into the group (r->G)
real(dp) :: greal(3)   ! Floating point representation of the group element components (g-vector)
integer  :: g(3)       ! Integer version of greal
integer  :: gIndx(n*size(pBas,2)) ! (ordinal) index of g-vector in the group
real(dp) :: Ainv(3,3) ! Inverse of the parent lattice vectors
logical err

stop "This routine is buggy!"
call matrix_inverse(pLV,Ainv,err)
if(err) stop "Coplanar lattice vectors in call to map_enumStr_to_real_space"
nD = size(pBas,2)
! Define the non-zero elements of the HNF matrix
a = HNF(1,1); b = HNF(2,1); c = HNF(2,2)
d = HNF(3,1); e = HNF(3,2); f = HNF(3,3)
! Compute the superlattice vectors 
sLV = matmul(pLV,HNF)

! Find the coordinates of the basis atoms
allocate(aBas(3,n*nD))
allocate(spin(n*nD))

! Map HNF points to real space
call matrix_inverse(real(HNF,dp),Sinv)

gTab2Latt2Cart = matmul(sLV,Sinv)

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
         aBas(:,ic) = matmul(gTab2Latt2Cart,(/z1,z2,z3/))+pBas(:,iD) 
      enddo
   enddo
enddo
enddo
if (ic /= n*nD) stop "ERROR: map_enumStr_to_real_space: Didn't find the correct # of basis atoms"
! Now map each position into the group so that the proper label can be applied
map2G = matmul(L,Ainv) ! r->G transformation
do iAt = 1,nD*n ! Loop over every atom
   greal = matmul(map2G,aBas(:,iAt)) ! Map position into the group
   g = nint(greal) ! Convert the g-vector from real to integer
   if(.not. equal(greal,g,eps)) stop "map2G didn't work in map_enumStr_to_real_space"
   ! Map the group element components (and d-vector index) to a single number
   ! This indexes the group member in the labeling
   gIndx(iAt) = (iAt-1)*n*nD+g(1)*S(2)*S(3)+g(2)*S(3)+g(3)+1
enddo

allocate(x(k))
x = 0.0
if (mod(k,2)==0) then
   do iAt = 1, n*nD
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
END MODULE enumeration_utilities
