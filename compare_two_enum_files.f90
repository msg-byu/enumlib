! This program is used to test different versions of the enum code---it compares the output of to
! make sure that the same structures are present in both files, even if the representations are not
! identical.
PROGRAM compare_2_enum_files
use num_types
use enumeration_utilities
implicit none

character(80) :: fname1,fname2,dummy
integer i1, i2, status

logical strFound 

if(iargc()/=2) stop "Need two arguments: name of each file to be compared"

call getarg(1,dummy)
read(dummy,'(a80)') fname1
call getarg(2,dummy)
read(dummy,'(a80)') fname2

open(11,file=fname1,status="old",iostat=status)
if(status/=0) stop "Filename does not exist"
open(12,file=fname2,status="old",iostat=status)
if(status/=0) stop "Filename does not exist"

i1 = 0
! Read in basic info from each file
do ! Loop over all of the structures in the first file
   i1 = i1 + 1
   i2 = 0
   do; i2 = i2 + 1 ! Loop over all the structures in the second file
      
   enddo
enddo

END PROGRAM compare_2_enum_files
