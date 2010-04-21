program compare_two_struct_enum
use num_types
use enumeration_types
use numerical_utilities
use enumeration_utilities
use io_utils
use derivative_structure_generator
implicit none
character(80) f1name, title, f2name, dummy
real(dp), dimension(3,3) :: pLV1, pLV2
integer,  dimension(3,3) :: L
real(dp), pointer :: dset1(:,:), dset2(:,:), sLVlist(:,:,:)
integer, pointer :: pLabel(:,:), HNF(:,:,:), eq(:), digit(:)
integer, pointer :: HNFin(:,:,:), label(:,:)
integer :: LatDim1, LatDim2, match, nD1, nD2, Nmin, Nmax, k, ioerr
integer :: strN, sizeN, n, pgOps, diag(3), a, b, c, d, e, f
real(dp) :: eps
logical full
character(maxLabLength)   :: labeling ! List, 0..k-1, of the atomic type at each site
type(RotPermList)         :: dRotList
type(RotPermList),pointer :: LattRotList(:)
type(opList), pointer     :: fixOp(:)


allocate(HNFin(3,3,1))
eps = 1e-4
print *, "Epsilon is currently hardwired at", eps
call getarg(1,dummy)
read(dummy,'(a80)') f1name
call getarg(2,dummy)
read(dummy,'(a80)') f2name
if(iargc()/=2) stop "Need two arguments: source and target struct_enum.out-type files"
call read_input(title,LatDim1,pLV1,nD1,dset1,k,eq,Nmin,Nmax,eps,full,pLabel,digit,f1name)
call read_input(title,LatDim2,pLV2,nD2,dset2,k,eq,Nmin,Nmax,eps,full,pLabel,digit,f2name)
if (LatDim1/=LatDim2) stop "Bulk/surf modes are not the same in the input files"
if (.not. equal(pLV2,pLV1,eps)) stop "Parent lattice vectors are not equivalent"
if (nD1/=nD2) stop "Number of d-vectors are not equivalent"
if (.not. equal(dset1,dset2,eps)) stop "D-vectors are not equivalent"
allocate(label(1,nD1),digit(nD1))
label = 1
digit = 1

open(10,file=f1name,status="old")
do ! read f1name file until the structure list begins
   read(10,*) title
   title = adjustl(title)
   if(title(1:5).eq."start")exit
enddo
print *,"Successfully opened first file and advanced to the structures list"
call get_dvector_permutations(pLV1,dset1,dRotList,LatDim1,eps)


do ! Read each structure from f1 and see if it is in the list of f2 structures
   read(10,*,iostat=ioerr) strN, sizeN, n, pgOps, diag, a,b,c,d,e,f, L, labeling
   if(ioerr/=0) exit
   HNFin = 0
   HNFin(1,1,1) = a; HNFin(2,1,1) = b; HNFin(2,2,1) = c;
   HNFin(3,1,1) = d; HNFin(3,2,1) = e; HNFin(3,3,1) = f;
   call remove_duplicate_lattices(HNFin,LatDim1,pLV1,dset1,dRotList,HNF,fixOp,&
                                  LattRotList,sLVlist,label,digit,eps)

!   call find_match_in_structenumout(f2name,pLV1,dset1,HNF,SNF,LatDim,pLabel,match,eps) 
   if(match/=0) then
      write(*,'("Structure #: ",i8," is a match to the structure in ",a80)') match, adjustl(f2name)
   else
      write(*,'("No match to the structure in ",a80)') adjustl(f2name)
      stop
   endif
   stop
enddo
endprogram compare_two_struct_enum

