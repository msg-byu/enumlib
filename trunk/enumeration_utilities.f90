MODULE enumeration_utilities
use num_types
use utilities_module
use symmetry_module
use vector_matrix_utilities
use numerical_utilities
use rational_mathematics
use enumeration_types
use derivative_structure_generator
use labeling_related
implicit none
private
public  map_enumStr_to_real_space, cartesian2direct, read_poscar, get_HNF_of_derivative_structure, &
        get_gspace_representation, find_match_in_structenumout

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
integer iAt, i, iD, nD
real(dp) :: greal(3)   ! Floating point representation of the group element components (g-vector)
integer  :: g(3)       ! Integer version of greal

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
         aBas(:,ic)=matmul(pLV,(/z1,z2,z3/))+pBas(:,iD)
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

!***************************************************************************************************
SUBROUTINE read_poscar(fname,title,sLV,aBas,aTyp)
character(80), intent(in) :: fname
character(80), intent(out):: title
real(dp), intent(out) :: sLV(3,3) ! parent and superlattice vector of the read-in structure
real(dp), pointer :: aBas(:,:) ! (out) the atomic basis vectors of the atoms (in Cartesian, not direct coords)
integer, pointer :: aTyp(:) ! (out) the atomic label (ie, type) for each atom in the atomic basis

integer iAt, iLV, status, ic, nTyp, nAt, idx
real(dp) scale, sLVinv(3,3)
character(80) coords, dummy, d2 ! Direct or Cartesian; input_string
logical err
integer, allocatable :: aList(:)

open(16,file=fname,iostat=status)
if(status/=0)then; write(*,'("The file specified for read_poscar in enumeration_utilities does not exist")')
   write(*,'("Filename: ",a80)') fname; stop; endif
open(17,file="readcheck_poscar_check.out")

read(16,'(a80)') title
write(17,'(a80)') title
read(16,*) scale
write(17,'("Scale factor: ",f8.5,/)') scale
read(16,*) (sLV(:,iLV),iLV=1,3)
write(17,'("Vectors (as columns)",/,3(3(f7.3,1x),/))') (sLV(:,iLV),iLV=1,3)
sLV = sLV*scale
write(17,'("Vectors (scaled)",/,3(3(f7.3,1x),/))') (sLV(:,iLV),iLV=1,3)
call matrix_inverse(sLV,sLVinv,err)
if (err) stop "Co-planar vectors in call to read_poscar"
read(16,'(a80)') dummy ! list of # of each atom type

print *, "Would be nice to put this in utilities instead or re-inventing it again elsewhere"
! Need to parse this line into a bunch of numbers and store them in an array
dummy=trim(adjustl(dummy)) ! Get rid of leading and trailing blanks
d2 = dummy
! Now find out how many number there are in the dummy line
ic = 0
do
   ic = ic + 1
   read(d2,*) status
   d2=d2(index(d2," "):80)
   d2 = adjustl(d2)
   if(len_trim(d2)==0) exit
   if(ic>25) stop "too many atoms in POSCAR: read_poscar in enumeration_utilities" 
enddo
nTyp = ic
allocate(aList(nTyp))
write(17,'(" Number of each kind of atom:")')
do ic = 1, nTyp
   read(dummy,*) aList(ic)
   dummy = dummy(index(dummy," "):80)
   write(17,'(i3)',advance="no") aList(ic)
enddo
write(17,*)
nAt = sum(aList)
allocate(aBas(3,nAt),aTyp(nAt))
write(17,'("Number of atoms in structure: ",i3)') nAt

read(16,*) coords
call ucase(coords)
write(17,'("Coordinates (direct or Cartesian): ",a1)') coords
! Read in all the atoms in the POSCAR. Convert to Cartesian coordinates, if necessary
do iAt = 1, nAt
   read(16,*) aBas(:,iAt)
   write(17,'("Atom #: ",i3," position: ",3(f7.3,1x))') iAt,aBas(:,iAt) 
   if(coords=="D") then ! convert to Cartesian
      aBas(:,iAt) = matmul(sLV,aBas(:,iAt))
      write(17,'("Atom #: ",i3," position (Cart): ",3(f7.3,1x))') iAt,aBas(:,iAt)
   endif
enddo
close(16)
! Set the type for each atom
idx = 1
do iAt = 1, nAt
   aTyp(iAt) = idx-1
   write(17,'("#:",i3,"   idx:",i3)') iAt, idx
   if(iAt == sum(aList(1:idx))) idx = idx + 1
enddo
close(17)
END SUBROUTINE read_poscar

!***************************************************************************************************
! This subroutine takes a structure defined in real space (atomic basis vectors in Cartesian 
! coordinates) and checks whether or not it is a derivative structure of a parent lattice. If it is,
! the parent lattice, d-set, and HNF/SNF are returned
SUBROUTINE get_HNF_of_derivative_structure(sLV,aBas,aTyp,pLV,dset,HNF,SNF,L,eps)
real(dp), intent(in) :: sLV(3,3), aBas(:,:) ! Input lattice vectors and atoms
integer, intent(in) :: aTyp(:) ! Type of each atom
real(dp), intent(out):: pLV(:,:) ! lattice vectors of the parent lattice
real(dp), pointer :: dset(:,:) ! (out)  d-set of the parent lattice
real(dp), intent(in):: eps
integer, pointer :: HNF(:,:,:) ! (out) List of (rotationally-) equivalent HNFs
integer, intent(out) :: SNF(3,3), L(3,3) ! L is the left transform for SNF

integer nAt, iAt, nOp, iOp, row(6), iuq, nuq, nMaxList, EntMax, nMatch, iEntry, i
integer, pointer :: aTypTemp(:)
real(dp), pointer :: aBasTemp(:,:), sgrot(:,:,:), sgshift(:,:)
real(dp) :: pLVinv(3,3)
integer,dimension(3,3) :: T, S, R, tempH, newH
logical err, unique
integer, allocatable :: trow(:,:), vs(:), idx(:)


nAt = size(aBas,2)
if(nAt/=size(aTyp)) stop "Input to get_HNF_of_derivative_structure is inconsistent"

open(17,file="debug_get_HNF.out")
write(17,'(3/,"<<< Finding the g-space representation >>>",2/)')

do iAt = 1, nAt
   write(17,'("Atom #: ",i2," position: ",3(f7.3,1x))') iAt, aBas(:,iAt)
enddo

allocate(aTypTemp(nAt),aBasTemp(3,nAt))
aTypTemp = aTyp; aBasTemp = aBas
pLV = sLV ! Use pLV as a temporary
call make_primitive(pLV,aTypTemp,aBasTemp,.false.,eps)
if(nAt/=size(aTypTemp)) then;
   write(17,'(/,"ERROR: The input structure wasn''t primitive")')
   write(17,'("atom type: ",20(i2,1x))') aTypTemp(:)
   write(17,'("number of atoms: ",i2,5x," size of aTyp:",i2)') nAt, size(aTypTemp)
   stop "ERROR: input structure for get_HNF_of_derivative_structure was not primitive";endif
aTypTemp = 1; pLV = sLV
call make_primitive(pLV,aTypTemp,aBasTemp,.false.,eps)
write(17,'("Parent lattice: ",/,3(3(f7.3,1x),/))') (pLV(:,iAt),iAt=1,3)
allocate(dset(3,size(aTypTemp))); dset = aBasTemp
write(17,'("d-set: ",/,200(3(f7.3,1x),/))') (dset(:,iAt),iAt=1,size(dset,2))

call matrix_inverse(pLV,pLVinv,err)
if(err) stop "Coplanar vectors in get_HNF_of_derivative_structure"
write(17,'("Size of supercell: ",i3)') abs(nint(determinant(sLV)/determinant(pLV)))

if(.not. equal(determinant(sLV)/determinant(pLV),nint(determinant(sLV)/determinant(pLV)),eps)) stop &
     "The superlattice in not an integer multiple of the primitive"
if(.not. equal(matmul(pLVinv,sLV),nint(matmul(pLVinv,sLV)),eps)) stop &
     "ERROR: HNF was non-integer"
S = nint(matmul(pLVinv,sLV))

call get_spaceGroup(pLV,aTypTemp,dset,sgrot,sgshift,.false.,eps)
nOp = size(sgrot,3)
allocate(trow(nOp,6),vs(nOp),idx(nOp))
idx = (/(i,i=1,nOp)/)
nuq = 0
do iOp = 1, nOp
   ! Under the current rotation, what HNF does the current HNF turn into?
   tempH = nint(matmul(pLVinv,matmul(sgrot(:,:,iOp),matmul(pLV,S))))
   call HermiteNormalForm(tempH,newH,R)
   unique = .true.
   row = (/newH(1,1),newH(2,1),newH(2,2),newH(3,1),newH(3,2),newH(3,3)/)
   do iuq = 1, nuq
      if (all(row == trow(iuq,:))) then
         unique = .false.
         exit
      endif
   enddo
   if (unique) then
      nuq = nuq + 1
      trow(nuq,:) = row
   endif
enddo
write(17,'("Unique HNF entries after rotations:")')
do i = 1,nuq
   write(17,'("row #:",i3,3x,6(i2,1x))') i,trow(i,:)
enddo
allocate(HNF(3,3,nuq))
HNF = 0
do iuq = 1,nuq
   HNF(1,1,iuq) = trow(iuq,1); HNF(2,1,iuq) = trow(iuq,2); HNF(2,2,iuq) = trow(iuq,3)
   HNF(3,1,iuq) = trow(iuq,4); HNF(3,2,iuq) = trow(iuq,5); HNF(3,3,iuq) = trow(iuq,6)
   write(17,'("HNF #: ",i3,/,3(3(i2,1x),/))') iuq,(HNF(iAt,:,iuq),iAt=1,3)
enddo
!!! All of the unique HNFs generated by the symmetries of the parent lattice are equivalent. So
!!! although the each entry in the list of surviving HNFs is different, they represent the same
!!! (super)lattice. Given that we can have multiple representations (HNFs) for the same lattice, we
!!! need a way to order the representations, so that we can use just one during the comparisons. We'll
!!! do this by selecting the HNF with the largest entries (comparing right to left). We pick this HNF
!!! because it's the one that would appear in the list of HNFs generated by enumeration code in
!!! enumlib (so of the list, we know it's the one that should appear in the struct_enum.out file).
!!! Loop over the each of the 6 entries in the HNFs. Find those that are maximum and copy those rows
!!! to the top of the list. Quit when you know that the one at the top of the list is the one that
!!! will appear first in the struct_enum.out list. 
! Actually, the order that HNFs are elimintated in the enum code isn't what I thought it was
!  initially. For the sake of a robust code, better to pass out the entire list of HNFs that are
!  unique and check them all. We can rethink this if efficiency becomes an issue.

!nMaxList = nuq
!do iEntry = 6,1,-1
!   EntMax = maxval(trow(1:nMaxList,iEntry))
!   nMatch = count(trow(1:nMaxList,iEntry)==EntMax)
!   nMaxList = nMaxList - nMatch
!   vs(1:nMatch) = pack(idx(1:nuq),trow(1:nuq,iEntry)==EntMax)
!   trow(1:nMatch,:) = trow(vs(1:nMatch),:)
!   if (nMaxList==1) exit
!enddo
!write(17,'("Winning HNF  ",6i3)') trow(1,:)
!HNF = 0
!HNF(1,1) = trow(1,1); HNF(2,1) = trow(1,2); HNF(2,2) = trow(1,3)
!HNF(3,1) = trow(1,4); HNF(3,2) = trow(1,5); HNF(3,3) = trow(1,6)


call SmithNormalForm(S,L,SNF,R)
write(17,'("Integer transform of superlattice: ",/,3(3(i2,1x),/))') (S(iAt,:),iAt=1,3)
!write(17,'("HNF: ",/,3(3(i2,1x),/))') (HNF(iAt,:),iAt=1,3)
write(17,'("SNF: ",/,3(3(i2,1x),/))') (SNF(iAt,:),iAt=1,3)
close(17)
END SUBROUTINE get_HNF_of_derivative_structure


!***************************************************************************************************
! Map a list of real space atomic basis vectors and their labels into g-space and extracts the
! labeling 
SUBROUTINE find_labeling_from_atom_basis(L,A,aBas,aTyp,SNF,labeling,eps)
integer, dimension(3,3), intent(in) :: L
real(dp), dimension(3,3), intent(in) :: A
real(dp), dimension(:,:), pointer :: aBas ! (in)
integer, dimension(:), intent(in) :: aTyp
integer, dimension(3,3), intent(in) :: SNF
integer, pointer :: labeling(:) ! out
real(dp), intent(in) :: eps

integer, dimension(3,size(aBas,2)) :: g
integer, pointer :: p(:,:)
integer :: diag(3), perm(size(aBas,2))
integer :: i, iAt, nAt
real(dp) :: Ainv(3,3)
logical err

open(18,file="debug_map_vectors.out")
write(18,'(3/,"Takes a list of vectors, maps them into the group, determines the labeling",2/)')

call matrix_inverse(A,Ainv,err)
if(err) stop "Coplanar vectors in find_labeling_from_atom_basis"

diag(1) = SNF(1,1); diag(2) = SNF(2,2); diag(3) = SNF(3,3) 
nAt = size(aBas,2)

write(18,'("Diagonal entries of HNF: ",3(i3,1x))') diag
write(18,'("Left transform:",/,3(3i3,1x,/),/)') transpose(L)
write(18,'("parent lattice vectors (columns):",/,3(3(f7.3,1x),/),/)') transpose(A)
write(18,'("parLatt inverse (columns):",/,3(3(f7.3,1x),/),/)') transpose(Ainv)
do iAt = 1, nAt
   write(18,'("Atom #: ",i3,"   position: ",3(f7.3,1x))') iAt,aBas(:,iAt)
enddo

g = 0
do iAt = 1, nAt ! Map each real space vector (atom position) into the group
   g(:,iAt)=matmul(matmul(L,Ainv),aBas(:,iAt))
   g(:,iAt) = modulo(g(:,iAt),diag) ! Move each vector into first cell in g-space
enddo

! Print out g-space representation of each vector
write(18,'("group list:")')
do i = 1, 3
   write(18,'(200(i2,1x))') g(i,:)
enddo

! Find the original (unpermuted) group
call make_member_list(diag,p)
write(18,'("original member list:")')
do i = 1, 3
   write(18,'(200(i2,1x))') p(i,:)
enddo

! Find the permutation of the original group that gives the input atom types
call find_permutation_of_group(p,g,perm)
allocate(labeling(nAt))
labeling = aTyp(perm)
write(18,'("Input atom labels: ",200(i2,1x))') aTyp
write(18,'("Group order:       ",200(i2,1x))') perm
write(18,'("Labeling:          ",200(i2,1x))') labeling

ENDSUBROUTINE find_labeling_from_atom_basis

!***************************************************************************************************
! This routine takes a derivative structure and finds its g-space representation. The routine
! returns all of the equivalent labelings of the superlattice so that these can be used to compare
! against possible matches in struct_enum.out. 
SUBROUTINE get_gspace_representation(pLV,dset,sLV,aBas,aTyp,HNF,LatDim,pLabel,eps)
real(dp), intent(in) :: pLV(3,3) ! Parent lattice vectors
real(dp), pointer ::  dset(:,:)  ! (in) "atomic basis" vectors of the multilattice
real(dp), intent(in) :: sLV(3,3) ! Superlattice vectors
real(dp), pointer :: aBas(:,:)   ! Atom basis vectors
integer, intent(in) :: aTyp(:)   ! Labels for atoms in the basis of the superlattice
integer, intent(in) :: HNF(:,:,:)! HNFs to generate the superlattice from pLV
real(dp), intent(in) :: eps
integer, intent(in) :: LatDim
integer, pointer :: pLabel(:,:)  ! The list of permuted labels (all equivalent)
!logical err 

integer, pointer :: aTypTemp(:), SNFlabel(:)
real(dp), pointer :: aBasTemp(:,:), sLVlist(:,:,:)
real(dp) :: pLVtemp(3,3),sLVtemp(3,3)
integer iAt, nAt, iOp, j, ip, iuq, nuq, nP, SNF(3,3)
integer,pointer,dimension(:,:,:) :: HNFin, HNFout, L, R, SNFlist, uqSNF
type(RotPermList) :: dRotList ! This is a list of permutations for the d-set 
                              ! (needed as in put for several routines)
type(opList), pointer :: fixOp(:) ! List of symops that fix the superlattice
type(RotPermList), pointer :: LattRotList(:) ! List of rotation perms for the superlattice
integer, pointer :: labeling(:), tempLabeling(:,:)
logical unique

open(17,file="debug_gspace_rep.out")
nAt = size(aBas,2)
if(nAt/=size(aTyp)) stop "Input to get_gspace_representation is inconsistent: nAt"

allocate(HNFin(3,3,1))
HNFin(:,:,1) = HNF(:,:,1)

allocate(aTypTemp(nAt),aBasTemp(3,nAt))
aTypTemp = aTyp; aBasTemp = aBas
sLVtemp = sLV 
call make_primitive(sLVtemp,aTypTemp,aBasTemp,.false.,eps)
if(nAt/=size(aTypTemp)) then;
   write(17,'(/,"ERROR: The input structure wasn''t primitive")')
   write(17,'("atom type: ",20(i2,1x))') aTypTemp(:)
   write(17,'("number of atoms: ",i2,5x," size of aTyp:",i2)') nAt, size(aTypTemp)
   stop "ERROR: input structure for get_gspace_representation was not primitive";endif
aTypTemp = 1; pLVtemp = sLV
call make_primitive(pLVtemp,aTypTemp,aBasTemp,.false.,eps)
write(17,'("Parent lattice: ",/,3(3(f7.3,1x),/))') (pLVtemp(:,iAt),iAt=1,3)
if (.not. equal(pLV,pLVtemp,eps)) stop "Input for get_gspace_representation is inconsistent"
write(17,'("Size of supercell: ",i3)') abs(nint(determinant(sLV)/determinant(pLV)))
write(17,'("d-set: ",/,200(3(f7.3,1x),/))') (dset(:,iAt),iAt=1,size(dset,2))
!** Calls to enumlib routines **
  ! This call generates a list of permutations for the d-set under symmetry operations
  ! of the parent lattice (need this for d-g table permutations)
call get_dvector_permutations(pLV,dset,dRotList,LatDim,eps)
write(17,'("d-set permutations: ",200(i3,1x))') dRotList%RotIndx(:)
  ! This call returns a list of operations that fix the superlattice. The routine expects a *list* of
  ! HNF matrices, but here we only need to pass in one because every one in the list is
  ! rotationally-equivalent. So the input and output lists are only one element
  ! long (in the last index). E.g., HNFin is a 3x3x1 array of integers (original lattice HNF)
call remove_duplicate_lattices(HNFin,LatDim,pLV,dset,dRotList,HNFout,fixOp,LattRotList,sLVlist,eps)
write(17,'("Number of symmetry operations that fix the superlattice: ",i3,/)') size(fixOp(1)%rot,3)
write(17,'(200(3(3f7.3,1x,/),"shift:",3(f7.3,1x),//))') &
     ((fixOp(1)%rot(j,:,iOp),j=1,3),fixOp(1)%shift(:,iOp),iOp=1,size(fixOp(1)%rot,3))
  ! This routine gets the SNF form of the HNF, returns a list of permutations effected by the
  ! rotation symmetries,
call get_SNF(HNFout,L,SNFlist,R,LattRotList,uqSNF,SNFlabel,fixOp)
write(17,'("Left transform:",/,3(3i2,1x,/))') (L(:,j,1),j=1,3)
write(17,'("Smith Normal form:",/,3(3i2,1x,/))') (SNFlist(:,j,1),j=1,3)
write(17,'("Permutations:",/,8(24(i3,1x),/))') LattRotList(1)%RotIndx(:)
SNF = SNFlist(:,:,1)

! This loads up the "perm" element of LattRotList
call get_rotation_perms_lists(pLV,HNFout,L,SNFlist,fixOp,LattRotList,dRotList,eps)
write(17,'("Permutations:")')
nP = size(LattRotList(1)%perm,1)
do ip = 1, nP
   write(17,'("Perm #",i3,":",1x,200(i2,1x))') ip,LattRotList(1)%perm(ip,:)
enddo
call find_labeling_from_atom_basis(L(:,:,1),pLV,aBas,aTyp,SNF,labeling,eps)
write(17,'("Input atom labels: ",200(i2,1x))') aTyp
write(17,'("Labeling:          ",200(i2,1x))') labeling

! Now that we know the labeling for the input atomic basis, find all possible permutations of that
! labeling. These can be compared to the entries in struct_enum.out
iuq = 0; nuq = 0; allocate(tempLabeling(nP,nAt)); tempLabeling = 0
do ip = 1, nP
   unique = .true.
   do iuq = 1, nuq
      if(all(labeling(LattRotList(1)%perm(ip,:))==tempLabeling(iuq,:))) then
         unique = .false.
         exit
      endif
   enddo
   if (unique) then
      nuq = nuq + 1
      tempLabeling(nuq,:) = labeling(LattRotList(1)%perm(ip,:))
   endif
enddo
write(17,'(/,"Number of unique labelings: ",i5)') nuq
allocate(pLabel(nuq,nAt))
pLabel = tempLabeling(1:nuq,:)
do iuq = 1, nuq
   write(17,'("uq Labeling # :",i3,5x,"labeling:",1x,200(i1,1x))') iuq,pLabel(iuq,:)
enddo
close(17)
END SUBROUTINE get_gspace_representation

!***************************************************************************************************
! Reads in structure info from a struct_enum.out-type file and compares to the g-space
! representation of a structure.
SUBROUTINE find_match_in_structenumout(fname,pLV,dset,sLV,aBas,aTyp,HNFin,SNF,LatDim,pLabel,match,eps) 
character(80), intent(in) :: fname
real(dp), intent(in), dimension(3,3) :: pLV, sLV
integer, intent(in), dimension(3,3) :: SNF
real(dp), pointer :: dset(:,:), aBas(:,:) ! (in)
integer, pointer :: aTyp(:), pLabel(:,:) ! (in, out) atom types, permuted labelings
integer, pointer :: HNFin(:,:,:) ! (in) List of HNFs (all equivalent)
integer, intent(in) :: LatDim
integer, intent(out) :: match
real(dp) :: eps

character(80) title, bulksurf, dummy,labeling
real(dp) :: p(3,3), Ainv(3,3)
real(dp), allocatable :: dvec(:,:)
integer, allocatable :: ilabeling(:)
integer iStr, volume, HNF(3,3), iLab, nLab, ioerr, i, nD, k, strN, hnfN, sizeN
integer a, b, c, d, e, f, pgOps, diag(3), L(3,3), n, iHNF, nHNF
logical foundLab, foundHNF

open(13,file="readcheck_find_match_in.out")
write(13,'(5(/),"This file echoes the input for the find_match_in_structenumout routine and")')
write(13,'("prints some of the computed quantities as well. It is intended to")')
write(13,'("be useful for debugging and fixing bad input")')
write(13,'("***************************************************")')

volume = determinant(SNF)
!call read_nth_line_from_enumlist()
open(11,file=fname,status='old',iostat=ioerr)
if(ioerr/=0)then; write(*,'("Input file doesn''t exist:",a80)') trim(fname);stop;endif
! Read in the title from the struct_enum.out file
read(11,'(a80)') title; title = trim(title)
write(13,'("Title: ",a80)') title

! Read in surf/bulk mode marker
read(11,'(a1)') bulksurf
write(13,'("bulk/surf: ",a1)') bulksurf
call ucase(bulksurf)
if((LatDim==2 .and. bulksurf=="B") .or. (LatDim==3 .and. bulksurf=="S")) then 
   print *,"The bulk/surf parameter in "//trim(fname)//" not consistent"
   write(*,'("Bulk/surf parameter is: ",a1)') bulksurf
   write(*,'("LatDim input in call is: ",i1)') LatDim
   stop
endif

! Read in the parent lattice vectors
do i = 1,3; read(11,*) p(:,i); enddo
call matrix_inverse(p,Ainv)
write(13,'("Parent lattice:",/,3(3f7.3,1x,/))') (p(:,i),i=1,3) 
write(13,'("Parent lattice inverse matrix:",/,3(3f7.3,1x,/))') (Ainv(:,i),i=1,3) 
if (.not. equal(determinant(pLV),determinant(p),eps)) then
   write(13,'("Volume of input structure parent: ",f7.3)') determinant(pLV)
   write(13,'("Volume of struct_enum parent: ",f7.3)') determinant(p)
   !stop "Volumes of parent lattices don't match"
endif
!if (.not. equal(matmul(Ainv,pLV),nint(matmul(Ainv,pLV)),eps)) stop "Parent lattices don't match"

! Read in the number of d-vectors, then read each one in
read(11,*) nD
if(nD/=size(dset,2)) stop "Numbers of d-vectors don't match"
allocate(dvec(3,nD),ilabeling(nD*volume))
do i = 1,nD; read(11,*) dvec(:,i); enddo
write(13,'("d-vectors:",/,40(3f8.4,1x,/))') (dvec(:,i),i=1,nD) 


! Read in the number of labels, i.e., binary, ternary, etc.
read(11,'(i2)') k
write(13,'("Number of labels (k): ",i2)') k 
read(11,*)
read(11,*) eps
write(13,'("Finite precision parameter (epsilon): ",g17.10)') eps 

do 
   read(11,*) dummy
   if(dummy=="start") exit
enddo
write(13,'("Found ""start"" label. Now reading in structures line-by-line ")')  

nHNF = size(HNFin,3)

! Read in the info for each structure in the file
write(13,'("Looking for a structure of size (volume factor): ",i3)') volume
write(13,'("strN, hnfN, sizeN, n, pgOps, diag, a,b,c,d,e,f, L, labeling")')  
iStr = 0; match = 0
nLab = size(pLabel,1)
do 
   iStr = iStr + 1
   read(11,*,iostat=ioerr) strN, hnfN, sizeN, n, pgOps, diag, a,b,c,d,e,f, L, labeling
   if(ioerr/=0) exit
   read(labeling,'(50i1)') ilabeling(1:n*nD)
   L = transpose(L) ! Listed with columns as the fast index in the struct_enum.out but reads
                    ! in with rows as the fast index. 
   write(13,'(i11,1x,i7,1x,i8,1x,i2,1x,i2,1x,3(i3,1x),1x,6(i3,1x),1x,1x,50i1)') &
        strN, hnfN, sizeN, n, pgOps, diag, a,b,c,d,e,f,ilabeling
   ! Define the full HNF matrix
   HNF = 0
   HNF(1,1) = a; HNF(2,1) = b; HNF(2,2) = c
   HNF(3,1) = d; HNF(3,2) = e; HNF(3,3) = f
   if(n/=volume) then
      write(13,'("Volume doesn''t match for str #:",i8," ---Skipping")')  strN
      cycle
   endif
   foundHNF = .false.
   do iHNF = 1, nHNF
      if(.not. all(HNF==HNFin(:,:,iHNF))) cycle
      foundHNF = .true.
   enddo
   if(.not. foundHNF) then
      write(13,'("HNF doesn''t match for str #:",i8," ---Skipping")') StrN
      cycle
   endif
   foundLab = .false.
   do iLab = 1, nLab
      if(all(ilabeling==pLabel(iLab,:))) then
         foundLab = .true.
         if (match/=0) stop "BUG! Found more than one match in struct_enum.out file"
         match = iStr
         exit
      endif
   enddo
   if(.not. foundLab) then
      write(13,'("Labeling didn''t match for str #:",i8)') strN
   else
      write(13,'("Structure number:",i8," is a match!")') strN
   endif
enddo
close(13)
END SUBROUTINE find_match_in_structenumout
END MODULE enumeration_utilities
