! This is not a pretty program. Not for general distribution. Don't release it into the wild.
! It has really been used. It works but is fragile.
program compare_two_struct_enum
use num_types
use enumeration_types
use numerical_utilities
use enumeration_utilities
use io_utils
use derivative_structure_generator
implicit none
character(800) f1name, title, dummy
character(len=:), allocatable :: f2name
real(dp), dimension(3,3) :: pLV1, pLV2
integer, dimension(3,3,1):: L, SNF
real(dp), pointer :: dset1(:,:), dset2(:,:), sLVlist(:,:,:)
integer, pointer :: pLabel(:,:), HNFout(:,:,:), eq(:), digit(:), dlabel(:,:)
integer, pointer :: HNFin(:,:,:), label(:,:)
integer :: LatDim1, LatDim2, match, nD1, nD2, Nmin, Nmax, k, ioerr
integer :: strN, sizeN, n, pgOps, diag(3), a, b, c, d, e, f
integer :: iuq, nuq, iP, nP, lc, iStr2, HNFtest(3,3), hdgenfact, labdgen, totdgen
integer :: strN2, hnfN2, sizeN2, n2, pgOps2, diag2(3), iL, nL, idx2
integer, allocatable :: ilabeling(:), ilabeling2(:)
real(dp) :: eps
logical full, HNFmatch, foundLab
character(maxLabLength)   :: labeling ! List, 0..k-1, of the atomic type at each site
type(RotPermList)         :: dRotList
type(RotPermList),pointer :: LattRotList(:)
type(opList), pointer     :: fixOp(:)
integer, pointer          :: crange(:,:)
logical                   :: conc_check
integer, pointer          :: degen_list(:)
integer, pointer          :: rotProdLab(:)

SNF = 0
allocate(HNFin(3,3,1))
eps = 1e-4
write(*,'("Epsilon is currently hardwired at ", g11.4)') eps
call getarg(1,dummy)
read(dummy,'(a800)') f1name
call getarg(2,dummy)
read(dummy,'(a800)') f2name
if(iargc()/=2) stop "Need two arguments: source and target struct_enum.out-type files"
!call read_struct_enum_out(title,LatDim1,pLV1,nD1,dset1,k,eq,Nmin,Nmax,eps,full,pLabel,digit,f1name,crange,conc_check)
!!GLWH Seems like the 'plabel' variable was being used for two different purposes. I've renamed this
!!one 'dlabel'
call read_struct_enum_out_oldstyle(title,LatDim1,pLV1,nD1,dset1,k,eq,Nmin,Nmax,eps,full,dLabel,digit,f1name,crange)
!call read_struct_enum_out(title,LatDim2,pLV2,nD2,dset2,k,eq,Nmin,Nmax,eps,full,pLabel,digit,f2name,crange,conc_check)
!write(*,'("Read file ",a80)') f1name
call read_struct_enum_out(title,LatDim2,pLV2,nD2,dset2,k,eq,Nmin,Nmax,eps,full,dLabel,digit,crange,f2name)
!!print *,"latdim1",latdim1,"latdim2"
!!if (LatDim1/=LatDim2) stop "Bulk/surf modes are not the same in the input files"
if (.not. equal(pLV2,pLV1,eps)) stop "Parent lattice vectors are not equivalent"
if (nD1/=nD2) stop "Number of d-vectors are not equivalent"
!print *,dset1-dset2,eps
!print *,abs(dset2-dset1)>eps
!print *,equal(dset1,dset2,eps)
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
! If there are multiple d-vectors in the parent lattice, we need to know what permutations (if any)
! are allowed. The permutations are needed to eliminate duplicate lattices.

call get_dvector_permutations(pLV1,dset1,dRotList,LatDim1,eps)

print *, "Be aware that HNFs are directly compared"
print *, "Rotationally equivalent HNFs are not considered"
print *, "Change this in the future"

do ! Read each structure from f1 and see if it is in the list of f2 structures
   !read(10,*,iostat=ioerr) strN, sizeN, hdgenfact, labdgen, totdgen, idx, n,  pgOps, diag, a,b,c,d,e,f, L(:,:,1), labeling
   read(10,*,iostat=ioerr) strN, sizeN, n,  pgOps, diag, a,b,c,d,e,f, L(:,:,1), labeling
   if(ioerr/=0) exit
   L(:,:,1) = transpose(L(:,:,1)) ! Written out column-wise but read in row-wise. So fix it by
   ! transposing
   SNF(1,1,1) = diag(1); SNF(2,2,1) = diag(2); SNF(3,3,1) = diag(3)

   HNFin = 0 ! Load up the HNF with the elements that were read in.
   HNFin(1,1,1) = a; HNFin(2,1,1) = b; HNFin(2,2,1) = c;
   HNFin(3,1,1) = d; HNFin(3,2,1) = e; HNFin(3,3,1) = f;
   !print *,"before remove"
   call remove_duplicate_lattices(HNFin,LatDim1,pLV1,dset1,dRotList,HNFout,fixOp,&
                                  LattRotList,sLVlist,degen_list,eps)
   open(17,file="debug_rotation_permutations.out")

   ! Get the list of label permutations
   !print *,"before got_rotperm"
   call get_rotation_perms_lists(pLV1,HNFout,L,SNF,fixOp,LattRotList,dRotList,eps)
   write(17,'("Rots Indx:",/,8(24(i3,1x),/))') LattRotList(1)%RotIndx(:)
   write(17,'("Permutation group (trans+rot):")')
   nP = size(LattRotList(1)%perm,1)
   do ip = 1, nP
      write(17,'("Perm #",i3,":",1x,200(i2,1x))') ip,LattRotList(1)%perm(ip,:)
   enddo

   ! Use the permutations effected by the rotations that fix the superlattice to generate labelings
   ! that are equivalent. The list of equivalent labelings will be used when we look for a match in the
   ! struct_enum file
   allocate(ilabeling(n*nD1),ilabeling2(n*nD1))
   read(labeling,'(500i1)') ilabeling
   !print *,"before find_eq"
   !print *,"ilabeling",ilabeling
   !print *,"size",size(ilabeling)
   !print *,"strN",strN
   call find_equivalent_labelings(ilabeling,LattRotList,pLabel,rotProdLab)
   !print *,"after find_eqv"
   nuq = size(pLabel,1)
   write(17,'(/,"Number of unique labelings: ",i5)') nuq
   do iuq = 1, nuq
      write(17,'("uq Labeling # :",i3,5x,"labeling:",1x,200(i1,1x))') iuq,pLabel(iuq,:)
   enddo
   close(17)
   match = 0
   ! Read in each structure from the second file and see if it matches the current structures from
   ! file 1.
   open(11,file=f2name,status="old")
   lc = 0 ! Count the number of lines
   do ! read f1name file until the structure list begins
      read(11,*) title
      title = adjustl(title)
      if(title(1:5).eq."start")exit
      lc = lc + 1
      if (lc > 100) stop "Didn't find the 'start' tag in the second file"
   enddo
   !print *,"Successfully opened second file and advanced to the structures list"
   open(13,file="debug_match_check.out")
   iStr2 = 0
   do
      iStr2 = iStr2 + 1
!      read(11,*,iostat=ioerr) strN2, hnfN2, sizeN2, n2, pgOps2, diag2, a,b,c,d,e,f, L, labeling
      read(11,*,iostat=ioerr) strN2, hnfN2, hdgenfact, labdgen, totdgen,  idx2, n2,  pgOps2, diag2, a,b,c,d,e,f, L(:,:,1), labeling
!      print *,"hnf",a,b,c,d,e,f,"lab2",labeling

      if(ioerr/=0) exit
      if (n2 < n) then ! Unit cells in this block are too small
         write(13,'("Volume is too small for structure #: ",i9," in file 2 (label # ",i9,")")') iStr2
         cycle 
      endif
      if (n2 > n) then  ! We've passed the point in the f2 file where the size of cells matches
         write(13,'("Volume is too big for structure #: ",i9," in file 2 (label # ", &
              & i9,")")') iStr2, strN2
         !exit
      else
         write(13,'("Volume matches for structure #: ",i9," in file 2 (label # ", &
              & i9,")")') iStr2, strN2
      endif
      read(labeling,'(500i1)') ilabeling2(1:n2*nD2)
!      print *,"ilab2",ilabeling2
      HNFtest = 0; 
      HNFtest = reshape((/a,b,d,0,c,e,0,0,f/),(/3,3/))
      !print *,HNFtest
      HNFmatch = .false.
      if(all(HNFtest==HNFin(:,:,1))) then ! HNFs match, next check the labeling
         write(13,'("HNF matches for structure #: ",i9," in file 2 (label # ", &
              & i9,")")') iStr2, strN2
         HNFmatch = .true.
      else
         write(13,'("HNF doesn''t match for structure #: ",i9," in file 2 (label # ", &
                       & i9,")")') iStr2, strN2
         cycle
      endif
      foundLab = .false.
      nL = size(pLabel,1)
      print *,"size p",shape(plabel)
      do iL = 1, nL
         write(*,'("ilabel",10(i1))') ilabeling2(:)
         write(*,'("plabel",10(i1))') plabel(iL,:)
         if(all(ilabeling2==pLabel(iL,:))) then
            foundLab = .true.
            if (match/=0) stop "BUG! Found more than one match in struct_enum.out file"
            match = iStr2
            exit
         endif
      enddo
!      stop "here"
      if(.not. foundLab) then
         write(13,'("Labeling didn''t match for str #:",i8)') strN
      else
         write(13,'("Structure number:",i8," is a match!")') strN
      endif
   enddo
   if(match==0) then
      write(13,'("Match for str #:",i8," in file 1 was not found in file 2")') strN
      write(*,'("Match for str #:",i8," in file 1 was not found in file 2")') strN
   else
      write(13,'("Structure number:",i8," is a match to # ",i9," in file 2.")') strN, match
      write(*,'("Structure number:",i8," is a match to # ",i9," in file 2.")') strN, match
   endif
   close(13)

   deallocate(ilabeling2)
   !print *,allocated(ilabeling2)
   deallocate(ilabeling)
!   if(strN==3) stop "here2"
   !read(*,*)
enddo
endprogram compare_two_struct_enum

