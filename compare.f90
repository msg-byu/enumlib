PROGRAM compare
use num_types
use numerical_utilities, only: equal
implicit none

character(80) :: file1, file2
character(40) :: lab1, lab2
integer status, k1, k2, n1(2), n2(2), dummy, nAt1, nAt2, pg1, pg2, SNF1(3), SNF2(3)
integer HNF1(6), HNF2(6), L1(9), L2(9), nAtCur, c1, c2
integer iSt, jSt, i
real(dp), dimension(3,3) :: A1, A2
real(dp) eps
logical match

eps = 1e-5_dp
read(*,*) file1, file2

! Open the files for reading
open(10,file=file1,status='old',iostat=status,action='read')
if(status/=0)then; write(*,'("Couldn''t find input file: ",a80)')adjustl(file1);stop;endif
open(11,file=file2,status='old',iostat=status,action='read')
if(status/=0)then; write(*,'("Couldn''t find input file: ",a80)')adjustl(file2);stop;endif

! Get the preamble stuff
read(10,*) 
read(11,*) 
do i = 1,3 ! Read in the lattice vectors
   read(10,*) A1(:,i)
   read(11,*) A2(:,i)
enddo
if(.not. equal(A1,A2,eps)) stop "Lattice vectors for the two files are incompatible"
read(10,'(i2)') k1
read(11,'(i2)') k2
if(k1/=k2) stop "Number of labels is not the same in both input files"
read(10,*) n1
read(11,*) n2
read(10,*); read(10,*)
read(11,*); read(11,*)
nAtCur = n1(1)
!allocate(lab1(n1(1)),lab2(n1(1)))

do ! Loop over each structure in file1
   read(10,*,iostat=status) iSt, dummy, nAt1, pg1, SNF1, HNF1, L1, lab1
   if (status/=0) exit ! End of file1 so exit the loop
   write(*,'("Reading str #: ",i10)') iSt
   !if (nAt1/=nAtCur) then; nAtCur = nAt1; deallocate(lab1,lab2)
   !   allocate(lab1(nAtCur), lab2(nAtCur))
   !   backspace(10)
   !   read(10,*) iSt, dummy, nAt1, pg1, SNF1, HNF1, L1, lab1
   !endif
   ! We now have a structure read from file1. Need to find a match in file2
   ! Back up file 2 until we are in the preceding SNF block
   do 
      backspace(11)
      read(11,*,iostat=status) jSt, dummy, nAt2, pg2, SNF2
      !write(*,*) SNF2,snf1
      if (status/=0) exit
      if (any(SNF1/=SNF2)) exit
      backspace(11)
   enddo
   ! Now pass through the file until we find a matching HNF
   match = .false.
   do 
      read(11,*,iostat=status) jSt, dummy, nAt2, pg2, SNF2, HNF2, L2, lab2
      write(*,'(i10,1x,i2,1x,i2,3(i3),1x,6(i3),2x,a20)') jSt, nAt2, pg2, SNF2, HNF2, adjustl(lab2)
      if (status/=0) exit
      if (nAt1/=nAt2) exit
      if (pg1/=pg2) cycle
      if (any(SNF1/=SNF2)) cycle
      if (any(HNF1/=HNF2)) cycle ! same must be true for L's. No need to check
      if (lab1/=lab2) then ! Need to check for rotation equivalence
         c1 = 0; c2 = 0
         do i = 1,len(lab1)
            if (lab1(i:i).eq.'0') c1 = c1 + 1
            if (lab2(i:i).eq.'0') c2 = c2 + 1
         enddo
         if (c1==c2) then
            write(*,'("Possible match: ",2i10,/,a20,/,a20)')&
                 iSt, jSt, adjustl(lab1), adjustl(lab2)
         endif
         cycle
      endif
      match = .true.
      exit
   enddo
   if (.not. match) then; write(*,'("No match for structure # ",i10," found in ",a80)')&
      iSt, adjustl(file2); stop
   endif
   write(*,'("Str. #: ",i10," matches ",i10)') iSt, jSt
   
enddo


close(10)
close(11)

ENDPROGRAM compare
