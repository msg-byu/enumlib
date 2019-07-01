! This is not a pretty program. Not for general distribution. Don't release it into the wild.
! It has really been used and it works...but is fragile.
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
implicit none
character(800)  title, dummy
character(len=:), allocatable :: f1name, f2name
real(dp), dimension(3,3) :: pLV1, pLV2, L2inv, T!, Tinv
integer, dimension(3,3,1):: L1, SNF, L2
real(dp), allocatable :: dset1(:,:), dset2(:,:), rot(:,:,:), shift(:,:)
real(dp), pointer :: sLVlist(:,:,:)
integer, pointer :: HNFout(:,:,:), eq(:), digit(:), HNFin(:,:,:)
integer, pointer, dimension(:,:) :: label, g1, g2, dlabel, pLabel, qLabel
integer :: LatDim1, LatDim2, match, nD1, nD2, Nmin, Nmax, k, ioerr
integer :: strN, sizeN, n, pgOps, diag(3), a, b, c, d, e, f, i, jl
integer :: iuq, nuq, iP, nP, lc, iStr2, HNFtest(3,3), hdgenfact, labdgen, totdgen
integer :: strN2, hnfN2, n2, pgOps2, diag2(3), iL, nL, idx, f1, f2
!integer :: sizeN2
integer, allocatable :: ilabeling(:), ilabeling2(:), atomType(:)
real(dp) :: eps
logical full, HNFmatch, foundLab, err
character(maxLabLength)   :: labeling ! List, 0..k-1, of the atomic type at each site
type(RotPermList)         :: dRotList1, dRotList2
type(RotPermList),pointer :: LattRotList1(:), LattRotList2(:)
type(opList), pointer     :: fixOp(:)
integer, pointer          :: crange(:,:)
!logical                   :: conc_check
integer, pointer          :: degen_list(:)
integer, pointer          :: rotProdLab(:)


SNF = 0
allocate(HNFin(3,3,1))
eps = 1e-4
write(*,'("Epsilon is currently hardwired at ", g11.4)') eps

! Read command line arguments for names of files to be compared
call getarg(1,dummy)
f1name = trim(dummy)
call getarg(2,dummy)
f2name = trim(dummy)
read(dummy,'(a800)') f2name
if(iargc()/=2) stop "Need two arguments: source and target struct_enum.out-type files"
!! This old "read" line was for a previous format of struct_enum.out, presumably.
!call read_struct_enum_out(title,LatDim1,pLV1,nD1,dset1,k,eq,Nmin,Nmax,eps,full,pLabel,digit,f1name,crange,conc_check)
!!GLWH Seems like the 'plabel' variable was being used for two different purposes. I've renamed this
!!one 'dlabel'
!call read_struct_enum_out_oldstyle(title,LatDim1,pLV1,nD1,dset1,k,eq,Nmin,Nmax,eps,full,dLabel,digit,f1name,crange)

print*,"Begin reading structure files"
call read_struct_enum_out(title,LatDim1,pLV1,nD1,dset1,k,eq,Nmin,Nmax,eps,full,dlabel,digit,cRange,f1name)
!call read_struct_enum_out(title,LatDim2,pLV2,nD2,dset2,k,eq,Nmin,Nmax,eps,full,pLabel,digit,f2name,crange,conc_check)
!write(*,'("Read file ",a80)') f1name
call read_struct_enum_out(title,LatDim2,pLV2,nD2,dset2,k,eq,Nmin,Nmax,eps,full,dLabel,digit,crange,f2name)
if (LatDim1/=LatDim2) stop "Bulk/surf modes are not the same in the input files"
print*,"Finished reading header info in both files"

! Basic checks for possible equivalence
if (.not. equal(pLV2,pLV1,eps)) stop "Parent lattice vectors are not equivalent"
if (nD1/=nD2) stop "Number of d-vectors are not equivalent"
if (.not. equal(dset1,dset2,eps)) stop "D-vectors are not equivalent"
allocate(label(1,nD1),digit(nD1))
label = 1
digit = 1

!"Fast forward" to structures in file1
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

open(13,file="debug_match_check.out")
do ! Read each structure from f1 and see if it is in the list of f2 structures
   read(f1,*,iostat=ioerr) strN, sizeN, hdgenfact, labdgen, totdgen, idx, n,  pgOps, diag, a,b,c,d,e,f, L1(:,:,1), labeling
   !read(f1,*,iostat=ioerr) strN, sizeN, n,  pgOps, diag, a,b,c,d,e,f, L(:,:,1), labeling
   if(ioerr/=0) exit
   write(*,'("file1 str: ",i5)') strN
   L1(:,:,1) = transpose(L1(:,:,1)) ! Written out column-wise but read in row-wise. So fix it by transposing
   SNF(1,1,1) = diag(1); SNF(2,2,1) = diag(2); SNF(3,3,1) = diag(3)

   HNFin = 0 ! Load up the HNF with the elements that were read in.
   HNFin(1,1,1) = a; HNFin(2,1,1) = b; HNFin(2,2,1) = c;
   HNFin(3,1,1) = d; HNFin(3,2,1) = e; HNFin(3,3,1) = f;

! Why are we doing this? To get the dRotList? We only have one HNF, so there are no duplicates to remove...
   call remove_duplicate_lattices(HNFin,LatDim1,pLV1,rot,shift,dset1,dRotList1,HNFout,fixOp,&
                                  LattRotList1,sLVlist,degen_list,eps)
   open(17,file="debug_rotation_permutations.out")

   ! Get the list of label permutations
   call get_rotation_perms_lists(pLV1,HNFout,L1,SNF,fixOp,LattRotList1,dRotList1,eps)
   write(17,'("Rots Indx:",/,8(24(i3,1x),/))') LattRotList1(1)%RotIndx(:)
   write(17,'("Permutation group (trans+rot):")')
   nP = size(LattRotList1(1)%perm,1)
   do ip = 1, nP
      write(17,'("Perm #",i3,":",1x,200(i2,1x))') ip,LattRotList1(1)%perm(ip,:)
   enddo

   ! Use the permutations effected by the rotations that fix the superlattice to generate labelings
   ! that are equivalent. The list of equivalent labelings will be used when we look for a match in the
   ! struct_enum file
   allocate(ilabeling(n*nD1),ilabeling2(n*nD1))
   read(labeling,'(500i1)') ilabeling
   !write(*,'("ilabeling:",/,30(i1),/)') ilabeling
   !print *,"size",size(ilabeling)
   !print *,"strN",strN
   call find_equivalent_labelings(ilabeling,LattRotList1,pLabel,rotProdLab)
   ! do i = 1,size(plabel,1)
   !    write(*,'("equiv_labelings <> ",100(i1))') pLabel(i,:)
   ! enddo
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
   f2 = 14;
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
   ! The trouble, subtle bug was that one of the routines from derStrGen was opening another file with unit number 11
   ! this was clashing with th unit 11 that as being opened here..
   compareloop: do
      write(*,'("Compare loop: file1 strN ",i6)') strN
      write(13,'("Compare loop: file1 strN ",i6)') strN
      iStr2 = iStr2 + 1
      !write(*,'("iStr2: ",i4)') iStr2
      read(f2,*,iostat=ioerr) strN2, hnfN2, hdgenfact, labdgen, totdgen, sizeN, n2, pgOps2, diag2, a,b,c,d,e,f, L2(:,:,1), labeling
!      read(f2,*,iostat=ioerr) strN2, hnfN2, hdgenfact, labdgen, totdgen,  idx2, n2,  pgOps2, diag2, a,b,c,d,e,f, L2(:,:,1), labeling
      write(*,'("read from file 2: strN2:",i5)') strN2
      if(ioerr/=0) exit
      read(labeling,'(500i1)') ilabeling2
      L2(:,:,1) = transpose(L2(:,:,1)) ! Written out column-wise but read in row-wise. So fix it by transposing
      if (any(L1/=L2)) then
         !make the qlabel list of permutations
         call get_dvector_permutations(pLV2,dset2,nD2,rot,shift,dRotList2,LatDim2,eps)
         call remove_duplicate_lattices(HNFin,LatDim2,pLV2,rot,shift,dset2,dRotList2,HNFout,fixOp,&
                                        LattRotList2,sLVlist,degen_list,eps)
         call get_rotation_perms_lists(pLV2,HNFout,L2,SNF,fixOp,LattRotList2,dRotList2,eps)
         call find_equivalent_labelings(ilabeling2,LattRotList2,qLabel,rotProdLab)
         open(17,file="debug_rotation_permutations.out",position="append")

         write(17,'("Rots Indx:",/,8(24(i3,1x),/))') LattRotList2(1)%RotIndx(:)
         write(17,'("Permutation group (trans+rot):")')
         nP = size(LattRotList2(1)%perm,1)
         do ip = 1, nP
            write(17,'("Perm #",i3,":",1x,200(i2,1x))') ip,LattRotList2(1)%perm(ip,:)
         enddo
         nuq = size(qLabel,1)
         write(17,'(/,"Number of unique labelings: ",i5)') nuq
         do iuq = 1, nuq
            write(17,'("uq Labeling # :",i3,5x,"labeling:",1x,200(i1,1x))') iuq,qLabel(iuq,:)
         enddo
         close(17)

      else
         if(associated(qlabel)) deallocate(qlabel)
         allocate(qlabel(1,n2*nD2))
         qlabel(1,:) = ilabeling2
      endif

      write(*,'("Comparing structures (f1 vs f2): ",2(i6,1x))') strN, strN2
      write(13,'("Comparing structures (f1 vs f2): ",2(i6,1x))') strN, strN2
      if (n2 < n) then ! Unit cells in this block are too small
         write(13,'("Volume is too small for structure #: ",i9," in file 2 (label # ",i9,")")') iStr2
         cycle
      endif
      if (n2 > n) then  ! We've passed the point in the f2 file where the size of cells matches
         write(13,'("Volume is too big for structure #: ",i9," in file 2 (label # ", &
              & i9,")")') iStr2, strN2
         exit
      else
         write(13,'("Volume matches for structure #: ",i9," in file 2 (label # ", &
              & i9,")")') iStr2, strN2
      endif
      !read(labeling,'(500i1)') ilabeling2(1:n2*nD2)
      HNFtest = 0;
      HNFtest = reshape((/a,b,d,0,c,e,0,0,f/),(/3,3/))
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
      do i = 0, k-1
         if (count(ilabeling2==i)/=count(ilabeling==i)) then
            write(13,'("Concentrations don''t match for color ",i1)') i
            write(*,'("Concentrations don''t match for color ",i1)') i
            cycle compareloop
         endif
      enddo
      write(*,'("Concentrations match for structures ",2(i6,1x))') strN, strN2
      write(13,'("Concentrations match for structures ",2(i6,1x))') strN, strN2
      if (any(L1/=L2)) then ! We must find the automorphism that connects the two different groups and adjust ilabeling2 accordingly
         print*,"hardwired change of labeling2 -- this is just temporary while testing/developing"
         write(*,'("labeling before: ",20(i1))') ilabeling2
         ilabeling2 = ilabeling2((/1,7,6,5,4,3,2,8,14,13,12,11,10,9/))
         write(*,'("labeling after:  ",20(i1))') ilabeling2
         if (associated(g1)) deallocate(g1)
         if (associated(g2)) deallocate(g2)
         !1) build the member list

         call make_member_list(diag,g1) ! This routine allocates g1
         do i = 1,3
            write(*,'("g1: ",20(i2,1x))') g1(i,:)
         enddo
         print*

         !2) find the group automorphism
         call matrix_inverse(real(L2(:,:,1),dp),L2inv,err)
         if (err) stop "L1inv is not defined"
         T = matmul(L1,L2inv(:,:,1))
         ! do i = 1,3
         !    write(*,'("L1:  ",3(i2,1x))') L1(i,:,1)
         ! enddo; print*
         ! do i = 1,3
         !    write(*,'("L2:  ",3(i2,1x))') L2(i,:,1)
         ! enddo
         do i = 1,3
            write(*,'(" T:  ",3(f7.3,1x))') T(i,:)
         enddo

         allocate(g2(3,size(g1,2)))
         g2 = -1
         g2 = nint(matmul(T,g1))
         !if (any(g2==-1)) stop "Bug in group member mapping"
         do i = 1,3
            write(*,'("g2: ",20(i2,1x))') g2(i,:)
         enddo
         print*
         do i = 1,size(qLabel,1)
            write(*,'("qL: ",20(i1))') qLabel(i,:)
         enddo
         print*

         !stop
         !3) change ilabeling2 into the equivalent ilabeling1
         !4) proceed as normal
      endif
         foundLab = .false.
         nL = size(pLabel,1)
         outer: do iL = 1, nL
            do jL = 1, size(qlabel,1)
               ilabeling2 = qLabel(jL,:)
               write(*,'("plabel: ",14(i1))') plabel(iL,:)
               write(*,'("qlabel: ",14(i1))') qlabel(jL,:)
               if(all(ilabeling2==pLabel(iL,:))) then
                  foundLab = .true.
                  print*,"MATCH"
                  if (match/=0) then
                     write(*,'("Additional match found: ",i6," == ",i6)') strN, strN2
                     ! write(*,'("In file 2, structure ",i6)') strN2
                     ! write(*,'("was equivalent to ",i6)') match
                     ! write(*,'("in file 1")')
                     !stop "BUG! Found more than one match in struct_enum.out file"
                  endif
      !            match = iStr2
                  match = strN2
                  exit outer
               endif
            enddo
         enddo outer
         if(.not. foundLab) then
            write(13,'("Labeling didn''t match for str #:",i8)') strN
         else
            write(13,'("Structure number:",i8," in file2 is a match to ",i8," in file1!")') strN2, strN
         endif
         if (match/=0) exit
   enddo compareloop
   close(f2)
   if(match==0) then
      write(13,'("Match for str #:",i8," in file 1 was not found in file 2")') strN
      write(*,'("Match for str #:",i8," in file 1 was not found in file 2")') strN
   else
      write(13,'("Structure number:",i8," is a match to # ",i9," in file 2.")') strN, match
      write(*,'("Structure number:",i8," is a match to # ",i9," in file 2.")') strN, match
   endif
   deallocate(ilabeling2)
   deallocate(ilabeling)
enddo
close(13)

endprogram compare_two_struct_enum
