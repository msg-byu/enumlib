! Compares two output files, usually "struct_enum.out" from enumlib (enum.x) to see if they contain
! equivalent lists.
!
! Started around 2008, used on and off, extensively revised in 2019 when the Smith Normal Form
! code in symlib (rational_mathematics) was updated to be more numerically stable (but
! unfortunately that affected many of the unit tests of enumlib). This program is used to compare
! output data from different versions of the code. The challenge is that data may often be
! equivalent but NOT identical.
!
! A further challenge is that the format of the input and output files change slowly over time and
! this code is fragile to that.
!
! This is not a pretty program. Nor is it efficient. Not for general distribution.
! Don't release it into the wild. It has really been used, a lot, and it works...but is fragile.
! GLWH 2019
program compare_two_struct_enum
use num_types
use symmetry
use enumeration_types
use numerical_utilities
use enumeration_utilities
use io_utils
use derivative_structure_generator
use vector_matrix_utilities
use labeling_related
use combinatorics
implicit none
character(800)  title, dummy
character(len=:), allocatable :: f1name, f2name
real(dp), dimension(3,3) :: pLV1, pLV2, L2inv, T!, L1inv!, Tinv
integer, dimension(3,3,1):: L1, SNF, L2
real(dp), allocatable :: dset1(:,:), dset2(:,:), rot(:,:,:), shift(:,:)
real(dp), pointer :: sLVlist(:,:,:)
integer, pointer :: HNFout(:,:,:), eq(:), digit(:), HNFin(:,:,:)
integer, pointer, dimension(:,:) :: label, g1=>null(), g2=>null(), dlabel, qLabel
integer :: LatDim1, LatDim2, match, nD1, nD2, Nmin, Nmax, k, ioerr
integer :: strN, sizeN, n, pgOps, diag(3), a, b, c, d, e, f, i, jl, j
integer :: iP, nP, lc, iStr2, HNFtest(3,3), hdgenfact, labdgen, totdgen
integer :: strN2, hnfN2, n2, pgOps2, diag2(3), idx, f1, f2, iLP, nLP
integer, allocatable :: ilabeling1(:), ilabeling2(:), atomType(:), automorphism(:)
integer, allocatable :: tempLab(:)
real(dp) :: eps
logical full, HNFmatch, foundLab, err, concMatch
character(maxLabLength)   :: labeling ! List, 0..k-1, of the atomic type at each site
type(RotPermList)         :: dRotList1, dRotList2
type(RotPermList),pointer :: LattRotList2(:)
type(opList), pointer     :: fixOp(:)
integer, pointer          :: crange(:,:), degen_list(:), rotProdLab(:), labPerms(:,:)

SNF = 0
allocate(HNFin(3,3,1))
eps = 1e-4
write(*,'("Epsilon is currently hardwired at ", g11.4)') eps

!
! Read command line arguments for names of files to be compared
!
call getarg(1,dummy)
f1name = trim(dummy)
call getarg(2,dummy)
f2name = trim(dummy)
read(dummy,'(a800)') f2name
if(iargc()/=2) stop "Need two arguments: source and target struct_enum.out-type files"
! This old "read" line was for a previous format of struct_enum.out, presumably.
!call read_struct_enum_out(title,LatDim1,pLV1,nD1,dset1,k,eq,Nmin,Nmax,eps,full,pLabel,digit,f1name,crange,conc_check)
!call read_struct_enum_out_oldstyle(title,LatDim1,pLV1,nD1,dset1,k,eq,Nmin,Nmax,eps,full,dLabel,digit,f1name,crange)

call read_struct_enum_out(title,LatDim1,pLV1,nD1,dset1,k,eq,Nmin,Nmax,eps,full,dlabel,digit,cRange,f1name)
!call read_struct_enum_out(title,LatDim2,pLV2,nD2,dset2,k,eq,Nmin,Nmax,eps,full,pLabel,digit,f2name,crange,conc_check)
call read_struct_enum_out(title,LatDim2,pLV2,nD2,dset2,k,eq,Nmin,Nmax,eps,full,dLabel,digit,crange,f2name)
if (LatDim1/=LatDim2) stop "Bulk/surf modes are not the same in the input files"
print*,"Finished reading header info in both files"

! Basic checks for possible equivalence
if (.not. equal(pLV2,pLV1,eps)) stop "Parent lattice vectors are not identical"
if (nD1/=nD2) stop "Number of d-vectors are not equivalent"
if (.not. equal(dset1,dset2,eps)) stop "D-vectors are not equivalent"
allocate(label(1,nD1),digit(nD1))
label = 1
digit = 1

!
! "Fast forward" to structures in file1
!
f1 = 10; open(f1,file=f1name,status="old")
do ! read f1name file until the structure list begins
   read(f1,*) title
   title = adjustl(title)
   if(title(1:5).eq."start")exit
enddo
print *,"Successfully opened first file and advanced to the structures list"


! If there are multiple d-vectors in the parent lattice, we need to know what permutations (if any)
! are allowed. The permutations are needed to eliminate duplicate lattices.
allocate(atomType(nD1)); atomType = 1
call get_spaceGroup(pLV1, atomType, dset1, rot, shift, .false., eps)
do i = 1, size(dset1,2)
write(*,'(3("dset1:   ",3(f10.6,1x)))') dset1(:,i)
enddo
call get_dvector_permutations(pLV1,dset1,nD1,rot,shift,dRotList1,LatDim1,eps)

print *, "Be aware that HNFs are directly compared"
print *, "Rotationally equivalent HNFs are not considered"
print *, "Change this in the future"

call get_permutations((/(i,i=0,k-1)/),labPerms)
nLP = size(labPerms,2)
write(*,'("label Perms: ",30i2)') (labPerms(i,:),i=1,size(labPerms,1))


!
! Loop over all structures in file1
!
open(13,file="debug_match_check.out")
f1loop: do ! Read each structure from f1 and see if it is in the list of f2 structures
   read(f1,*,iostat=ioerr) strN, sizeN, hdgenfact, labdgen, totdgen, idx, n,  pgOps, diag, a,b,c,d,e,f, L1(:,:,1), labeling
   !read(f1,*,iostat=ioerr) strN, sizeN, n,  pgOps, diag, a,b,c,d,e,f, L(:,:,1), labeling
   if(ioerr/=0) exit
   write(13,'("file1 str: ",i5)') strN
   L1(:,:,1) = transpose(L1(:,:,1)) ! Written out column-wise but read in row-wise. So fix it by transposing
   SNF(1,1,1) = diag(1); SNF(2,2,1) = diag(2); SNF(3,3,1) = diag(3)

   HNFin = 0 ! Load up the HNF with the elements that were read in.
   HNFin(1,1,1) = a; HNFin(2,1,1) = b; HNFin(2,2,1) = c;
   HNFin(3,1,1) = d; HNFin(3,2,1) = e; HNFin(3,3,1) = f;

   ! Use the permutations effected by the rotations that fix the superlattice to generate labelings
   ! that are equivalent. The list of equivalent labelings will be used when we look for a match in the struct_enum file
   if (allocated(tempLab)) deallocate(tempLab)
   allocate(ilabeling1(n*nD1),ilabeling2(n*nD1),tempLab(n*nD1))!,permutedLab(n*nD1)) 26 July, permutedLab not used
   read(labeling,'(500i1)') ilabeling1
   !write(*,'("ilabeling1:",/,30(i1),/)') ilabeling1
   match = 0
   ! Read in each structure from the second file and see if it matches the current structures from
   ! file 1.

!! END read from file1

   f2 = 14; ! Define a unit number for the second file.
   open(f2,file=f2name,status="old")
   lc = 0 ! Count the number of lines
   do ! read the second file until the structure list begins
      read(f2,*) title
      title = adjustl(title)
      if(title(1:5).eq."start")exit
      lc = lc + 1
      if (lc > 100) stop "Didn't find the 'start' tag in the second file"
   enddo
   iStr2 = 0

!! END reading of 'header' from file2
! You can get a subtle bug here, hard to chase down, if you have another module opening files with the same unit number. Make sure that the unit numbers for files opened by this program are unique.  GLWH 2019 (that's why f2 is set to 14 above)

   compareloop: do
      write(13,'("Compare loop: file1 strN ",i6)') strN
      iStr2 = iStr2 + 1
      read(f2,*,iostat=ioerr) strN2, hnfN2, hdgenfact, labdgen, totdgen, sizeN, n2, pgOps2, diag2, a,b,c,d,e,f, L2(:,:,1), labeling
      write(13,'("read from file 2: strN2: ",i5)') strN2
      if(ioerr/=0) exit
      read(labeling,'(500i1)') ilabeling2
      L2(:,:,1) = transpose(L2(:,:,1)) ! Written out column-wise but read in row-wise. So fix it by transposing

      write(*,'("Comparing structures (f1 vs f2): ",2(i6,1x))') strN, strN2
      write(13,'("Comparing structures (f1 vs f2): ",2(i6,1x))') strN, strN2
      if (n2 < n) then ! Unit cells in this block are too small
         write(13,'("Volume is too small for structure #: ",i9," in file 2 (label # ",i9,")")') iStr2
         cycle
      endif
      ! Note here that we are assuming the structures are listed in ascending order by size. If
      ! that is not the case, one cannot use this early exit to speed things up. As it stands now
      ! though, the enum code always writes in ascending size order.
      if (n2 > n) then  ! We've passed the point in the f2 file where the size of cells matches
         write(13,'("Volume is too big for structure #: ",i9," in file 2 (label # ", &
              !& i9,")")') iStr2, strN2
         exit
      else
         write(13,'("Volume matches for structure #: ",i9," in file 2 (label # ", &
              !& i9,")")') iStr2, strN2
      endif
      !read(labeling,'(500i1)') ilabeling2(1:n2*nD2)
      HNFtest = 0;
      HNFtest = reshape((/a,b,d,0,c,e,0,0,f/),(/3,3/))
      HNFmatch = .false.
      if(all(HNFtest==HNFin(:,:,1))) then ! HNFs match, next check the labeling
         write(13,'("HNF matches for structure #: ",i9," in file 2 (label # ", &
              !& i9,")")') iStr2, strN2
         HNFmatch = .true.
      else
         write(13,'("HNF doesn''t match for structure #: ",i9," in file 2 (label # ", &
               !        & i9,")")') iStr2, strN2
         cycle
      endif
      ! If ilabeling2 *or* one if its complements has the correct concentration, then we need to continue in the check.

      do iLP = 1, nLP ! Loop over all complements of the labels (stored in 'labPerms')
         do i = 0,k-1 ! Loop over each color to make the complement of the label
            where  (ilabeling2==i) tempLab = labPerms(iLP,i+1)
         enddo
         concMatch = .true.
         do i = 0, k-1
            if (count(tempLab==i)/=count(ilabeling1==i)) then
               concMatch = .false.
            endif
         enddo
         if (concMatch) exit
      enddo
      if (.not. concMatch) then
         write(13,'("Concentrations (and all complements) don''t match")')
         cycle compareloop
      endif
      write(13,'("Concentrations (or complements) match for structures:",2(i6,1x))') strN, strN2
! <<< The next 100 lines or so compares the labelings---everything else matches at this point >>>

   ! We do "remove_duplicate_lattices" even though there are no duplicates to remove because the routine also generates the vector subscripts for the label permutations and we need these to generate the list of equivalent labelings
   ! We could use dRotList1 is used in the "remove_..." and "get_r..." calls below because, at this point, we've established that dRotList1 and dRotList2 would be the same, and it would be more efficient to find dRotList1 outside of the compare loop
   call get_dvector_permutations(pLV2,dset2,nD2,rot,shift,dRotList2,LatDim2,eps) ! get dRotList2
   call remove_duplicate_lattices(HNFin,LatDim2,pLV2,rot,shift,dset2,dRotList2,HNFout,fixOp,&
                                  LattRotList2,sLVlist,degen_list,eps) ! get fixing operations
   !open(17,file="debug_rotation_permutations.out")
   ! Get the list of label permutations
   call get_rotation_perms_lists(pLV2,HNFout,L2,SNF,fixOp,LattRotList2,dRotList2,eps)
   !write(17,'("Rots Indx:",/,8(24(i3,1x),/))') LattRotList2(1)%RotIndx(:)
   !write(17,'("Permutation group (trans+rot):")')
   nP = size(LattRotList2(1)%perm,1)
   ! do ip = 1, nP
   !    write(17,'("Perm #",i3,":",1x,200(i2,1x))') ip,LattRotList2(1)%perm(ip,:)
   ! enddo
   call find_equivalent_labelings(ilabeling2,LattRotList2,qLabel,rotProdLab)
   if (any(L1/=L2)) then ! We must find the automorphism that connects the two different groups and adjust ilabeling2 accordingly
      print*,"L1/=L2"
      if (associated(g1)) deallocate(g1); if (associated(g2)) deallocate(g2)
      !1) build the member list
      call make_member_list(diag,g1) ! This routine allocates g1
      !2) find the group automorphism
      call matrix_inverse(real(L2(:,:,1),dp),L2inv,err)
      if (err) stop "Inverse of L_2 matrix is not defined"
      T = matmul(L1(:,:,1),L2inv)
      allocate(g2(3,size(g1,2)))
      !g2 = -1
      g2 = nint(matmul(T,g1))
      g2 = modulo(g2,spread(diag,2,n))

      ! Next, find the rearrangement between the two groups (the automorphism)
      if (allocated(automorphism)) deallocate(automorphism)
      allocate(automorphism(n2))
      call find_permutation_of_group(g1,g2,automorphism)
      !write(13,"('automorphism')")
      !write(13,"(50i1)") automorphism
      ! 3) Apply the automorphism to every permuted labeling (qLabels)
      do j = 1,size(qLabel,1)
         do i = 0, n*nD2-1, n ! Each block of 'n' needs a copy of the automorphism
            qLabel(j,1+i:n+i) = qLabel(j,automorphism+i)
         enddo
      enddo
   endif
   !4) proceed as normal

   foundLab = .false.
   ! We need to loop over all permutations (to get the full orbit) because we do not know
   ! which member of the orbit is the representative.
   ! We also need to check the complement of each labeling because any of the complements ! may appear in either list. GLWH 2019, see email chain, subject: "Monk?", with Rod Forcade starting May 16 2019, and subject "Eureka" starting Sep 5 2019
   do jL = 1, size(qlabel,1)
      !ilabeling2 = qLabel(jL,:)
      ! if(all(ilabeling2==ilabeling1)) then
      !    foundLab = .true.
      !    write(13,'("MATCH  Perm #: ",i3)') jl
      !    if (match/=0) then
      !       write(13,'("Additional match found: ",i6," == ",i6)') strN, strN2
      !        write(13,'("In file 2, structure ",i6)') strN2
      !        write(13,'("was equivalent to ",i6)') match
      !        write(13,'("in file 1")')
      !       stop "BUG! Found more than one match in struct_enum.out file"
      !    endif
      !    match = strN2
      !    exit
      ! endif
      do iLP = 1, nLP ! Loop over all complements of the labels (stored in 'labPerms').
         do i = 0,k-1 ! Loop over each color to make the complement of the label
            where  (qLabel(jL,:)==i) ilabeling2 = labPerms(iLP,i+1)
         enddo
         if(all(ilabeling1==ilabeling2)) then
            foundLab = .true.
            if (match/=0) then
               write(13,'("Additional match found: ",i6," == ",i6)') strN, strN2
               write(13,'("In file 2, structure ",i6)') strN2
               write(13,'("was also equivalent to ",i6)') match
               write(13,'("in file 1")')
               stop "BUG! Found more than one match in struct_enum.out file"
            endif
            match = strN2
            exit !iLP loop over complements
         endif
      enddo !iLP loop
      if (foundLab) exit ! loop over permutations. Some labelings are invariant under some permutations. So without the exit, we will find a duplicate
   enddo ! jL Loop over permutations
   if(.not. foundLab) then
      write(13,'("Labeling didn''t match for str #:",i8)') strN
   else
      write(13,'("Structure number:",i8," in file2 is a match to ",i8," in file1!")') strN2, strN
   endif
   !if (match/=0) exit
   enddo compareloop
   close(f2)
   if(match==0) then
      write(13,'("Match for str #:",i8," in file 1 was not found in file 2")') strN
      stop "Match not found!"
   else
      write(13,'("Structure number:",i8," is a match to # ",i9," in file 2.")') strN, match
   endif
   deallocate(ilabeling2)
   deallocate(ilabeling1)
enddo f1loop
close(13)
write(*,'("The labelings in the two lists have 1-1 correspondence.")')

endprogram compare_two_struct_enum
