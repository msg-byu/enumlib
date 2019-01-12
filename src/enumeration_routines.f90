MODULE enumeration_routines
use num_types
use symmetry, only: get_spaceGroup, rm_3d_operations
implicit none
private
public make_inactive_table, getSpaceGroup_activeSitesOnly
CONTAINS

!!<summary>Make a table of inactive sites. Col1: iD, Col2: label for site.</summary>
!!<comments>This routine makes a table so that we can distingish "inactive" sites from those that are actually needed for the enumeration. Inactive sites are those that can only take one kind of atom in the enumeration. If these are explicitly included in the enumeration, they cause an unnecessry combinatoric explosion. So we remove them for the enumeration and then put them back inside the enumerated cells ex post facto.
!!  equivalencies (1xN_sites)
!!                |-----|
!!            1   |  1  |
!!            2   |  2  |
!!            3   |  3  |
!!            4   |  4  |
!!            5   |  3  | (label on site 5 must match that on site 3)
!!                |-----| (no degree of freedom for site 5)
!!</comments>
SUBROUTINE make_inactive_table(k,nD, equivTbl,nDFull,labelFull,dFull,d,label,digit,digitFull,inactives)
integer, intent(in)  :: k ! Number of labels (i.e., binary, ternary, etc.)
integer, intent(out):: nD ! Number of "active" sites
integer, intent(in)  :: equivTbl(:) ! see above
integer, intent(in)  :: nDFull ! Total number of sites,
integer, allocatable :: digitFull(:)
real(dp), allocatable:: d(:,:), dFull(:,:)
integer, intent(inout):: labelFull(:,:) ! (NDfullxNlabels)
integer, allocatable :: label(:,:), digit(:)
integer, allocatable :: inactives(:,:)

integer              :: nInAct ! Number of inactive sites
integer              :: i, iD, jInactive ! counters

! Count the number of sites that are *not* equivalent to others
nD = count( (/(i,i=1,nDFull)/)==equivTbl)
! Compute the number of sites with only one label allowed (all
! -1's in the 'labelFull' table except one). These sites have no
! freedom and so we'll remove them from the enumeration.
nInAct = count((/(k-1==count(labelFull(:,i)==-1),i=1,nDFull)/))
nD = nD - nInAct

print*,"**labelFull**"
do i = 1, size(labelFull,2)
  write(*,'(6(i2,1x))') labelFull(:,i)
enddo

! Set up arrays for the subset of sites that have some
! configurational freedom.
allocate(d(3,nD),label(size(labelFull,1),nD), digit(nD))
allocate(inactives(nInAct,2))
print*,"Size of inactives",size(inactives)
nD = 0; jInactive = 0
do iD=1,nDFull
 if (iD==equivTbl(iD) .and. labelFull(2,iD)/=-1) then !
    !this dset member is a unique point and has more than one
    !label allowed (if it only had one allowed label, position
    !2 of the label table 'labelFull' would be -1).
    nD=nD+1
    d(:,nD)     = dFull(:,iD)
    label(:,nD) = labelFull(:,iD)
    digit(nD)   = digitFull(iD)
elseif (iD/=equivTbl(iD)) then
      ! this dset member is equivalent (concerning the
      ! enumeration!) to a different point.  => Force the label
      ! and digit arrays to be equal for equivalent sites
    labelFull(:,iD) = labelFull(:,equivTbl(iD))
    digitFull(iD)   = digitFull(equivTbl(iD))
 else
    print*,"Found inactive site"
 ! Sites with only one label will be left out of the dFull table but
 ! we need to keep a list of these sites so we can add them
 ! back in later when we print out the list.
    jInactive = jInactive + 1
    inactives(jInactive,:) = (/iD,labelFull(1,iD)/)
 endif
enddo
print*,"inactives>>>"
do i = 1, size(inactives,1)
  write(*,'(2(i1,1x))') inactives(i,:)
enddo

END SUBROUTINE make_inactive_table

SUBROUTINE getSpaceGroup_activeSitesOnly(pLV,nD,nDfull,dFull,inactives,latDim,rot,shift,eps)
real(dp), intent(in)  :: pLV(3,3)
integer, intent(in)   :: nD, nDFull
real(dp), allocatable :: dFull(:,:)
integer, allocatable  :: inactives(:,:)
integer, intent(in)  :: latDim ! Lat Dimensionality, i.e., surface or bulk
real(dp), allocatable:: rot(:,:,:), shift(:,:)
real(dp), intent(in):: eps

integer, allocatable :: aTyp(:)
integer :: iC, iD ! loop counters
integer :: status

allocate(aTyp(nDfull),STAT=status)
if(status/=0)stop "Allocation failed in get_dvector_permutations: aTyp"

aTyp = 1
aTyp(inactives(:,1)) = 2 ! inactive sites should never map to active, label distinctly
call get_spaceGroup(pLV,aTyp,dFull,rot,shift,.false.,eps)
if(latDim==2) call rm_3D_operations(pLV,rot,shift,eps)

END SUBROUTINE getSpaceGroup_activeSitesOnly


END MODULE
