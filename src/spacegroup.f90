! stand alone driver to read in lattice vectors and atomic basis
! and return the spacegroup rotations and shift.
PROGRAM spacegroup
use num_types
use symmetry
implicit none

real(dp), dimension(3,3) :: LV
real(dp), allocatable    :: d(:,:)
integer, allocatable     :: typ(:)
character(120)           :: LVfile, dFile
integer iD, nD
integer nArg

nArg = command_argument_count()
if (nArg)/= 2
stop "ERROR: requires two argments. File names for lattice vectors and for atomic basis vectors"

! Read the two file names for LV and d from the command line

! Cound the number of atoms, allocate atom positions and typ

! Call the spacegroup routine

! Write out the

END PROGRAM spacegroup
