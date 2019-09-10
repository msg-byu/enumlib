MODULE enumeration_routines
use num_types
use symmetry, only: get_spaceGroup, rm_3d_operations
use rational_mathematics, only: gcd
implicit none
private
public make_inactive_table, getSpaceGroup_activeSitesOnly, adjust_crange_for_inactive_sites
CONTAINS

  !!<summary>test</summary>
  !!<parameter name="concRanges" regular="true">Input for concentration ranges. This will be modified
  !! and passed back out</parameter>
  !!<parameter name="inactvTbl" regular="true">List of inactive sites (col1) and their labels (col2)
  !! in the full cell.</parameter>
  !!<parameter name="LabelsTbl" regular="true">Columns list the labels for each site. Rows iterate
  !!over colors</parameter>
  !!<parameter name="nDFull" regular="true">Number of sites in the parent cell, as input, no
  !! inactive sites removed </parameter>
  SUBROUTINE adjust_crange_for_inactive_sites(concRanges,inactvTbl,LabelsTbl,nDFull)
  integer :: concRanges(:,:), inactvTbl(:,:)
  integer, intent(in) :: LabelsTbl(:,:)
  integer, intent(in) :: nDFull ! number of d-vectors in parent cell
  integer iC, nC, iD ! Counter for colors, number of colors (often k elsewhere)
  integer nAct ! number of active sites
  integer nInAct ! number of sites that are inactive for color #iC
  integer denom
  !!<comments> For a write-up of the logic of this routine, see 'notes_cRangeAdjustment' in the
  !! support folder</comments>

!del  write(*,'("iD  lbl")')
!del  do iC = 1,size(inactvTbl,1)
!del    write(*,'(i3,1x,i3)') inactvTbl(iC,:)
!del  enddo
  nAct = nDFull - size(inactvTbl,1)
  nC = size(concRanges,1)
!del  write(*,'(i2," sites are active")') nAct
!del  write(*,'(i2," colors")')nC
!del  write(*,'("Size of Label table",2(i3,1x))') shape(LabelsTbl)
!del  write(*,'("Labels Table:")')
!del  do ic = 1, nC
!del    write(*,'(250(i2,1x))') LabelsTbl(ic,:)
!del    print*
!del  enddo

!del  print*,'Pack',pack((/(iD,iD=1,nDFull)/),LabelsTbl(2,:)==-1)
  do iC = 1, nC ! for each color, subtract the fraction of sites that are inactive, multiply by f
    ! If the current color appears only in columns of the labels table in columns that have only one
    ! color, then this color is "inactive"---it only appears on inactive sites; only by itself
    if (count(LabelsTbl(:,pack((/(iD,iD=1,nDFull)/),LabelsTbl(2,:)==-1)) == iC-1) == &
    & count(LabelsTbl == iC-1)) then
!del     write(*,'("** Color ",i2," is inactive")') iC-1
     concRanges(iC,1:2) = 0
     concRanges(iC,3) = 1 ! This isn't necessary but it's aesthetic
    cycle
   endif
!del    write(*,'("Color: ",i1)') iC-1
    nInAct = count(inactvTbl(:,2)==iC-1)
!del    write(*,'(i2," are inactive")')nInAct
    denom = concRanges(iC,3)
    concRanges(iC,1:2) = nDFull*concRanges(iC,1:2) - nInact*denom
    concRanges(iC,3) = concRanges(iC,3)*nDFull
    concRanges(iC,1:2) = concRanges(iC,1:2)*nDFull
    concRanges(iC,3) = concRanges(iC,3)*nAct
!del    write(*,'("concR: ",3(i3,","))') concRanges(iC,:)
    if (concRanges(iC,1) < 0) concRanges(iC,1)  = 0
    if (concRanges(iC,2) < 0) concRanges(iC,2)  = 0
    if (concRanges(iC,1) > concRanges(iC,3)) concRanges(iC,1)  = concRanges(iC,3)
    if (concRanges(iC,2) > concRanges(iC,3)) concRanges(iC,2)  = concRanges(iC,3)
    !  write(*,'("concR: ",3(i3,","))') concRanges(iC,:)
    concRanges(iC,:) = concRanges(iC,:)/gcd(concRanges(iC,:))
!del      write(*,'("concR: ",3(i3,","))') concRanges(iC,:)
  enddo


END SUBROUTINE adjust_crange_for_inactive_sites

!!<summary>Make a table of inactive sites. Col1: iD, Col2: label for site.</summary>
!!<comments>This routine makes a table so that we can distingish "inactive" sites from those that
!! are actually needed for the enumeration. Inactive sites are those that can only take one kind of
!! atom in the enumeration. If these are explicitly included in the enumeration, they cause an
!! unnecessry combinatoric explosion. So we remove them for the enumeration and then put them back
!! inside the enumerated cells ex post facto.
!! !
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
! nDFull is the total number of sites, nD is the number of non-equivalent sites
nD = count( (/(i,i=1,nDFull)/)==equivTbl)

! Compute the number of sites with only one label allowed (all
! -1's in the 'labelFull' table except the first one). These sites have no
! freedom and so we'll remove them from the enumeration.
! (Logic: for each row in the labels table, see how many of them have as many -1's as there
! are types of labels, minus 1)
nInAct = count((/(k-1==count(labelFull(:,i)==-1),i=1,nDFull)/))
! Now, nD will be the number of non-equivalent, non-inactive sites
nD = nD - nInAct

! Set up arrays for the subset of sites that have some
! configurational freedom.
allocate(d(3,nD),label(size(labelFull,1),nD), digit(nD))
allocate(inactives(nInAct,2))
! now nD will be used as a counter of the number of active sites
nD = 0; jInactive = 0
do iD=1,nDFull
  if (iD==equivTbl(iD) .and. labelFull(2,iD)/=-1) then !
    !this dset member is not an equivalent point and has more than one
    !label allowed (if it only had one allowed label, position
    !2 of the label table 'labelFull' would be -1).
    nD=nD+1
    d(:,nD)     = dFull(:,iD) ! Copy the position of the active site
    label(:,nD) = labelFull(:,iD) ! Copy the possible labels for the active site
    digit(nD)   = digitFull(iD) ! Copy the number of possible labels for this active site
  elseif (iD/=equivTbl(iD)) then
      ! this dset member is equivalent (concerning the
      ! enumeration!) to a different point.  => Force the label
      ! and digit arrays to be equal for equivalent sites
    labelFull(:,iD) = labelFull(:,equivTbl(iD))
    digitFull(iD)   = digitFull(equivTbl(iD))
  else
  ! Sites with only one label will be left out of the dFull table but
  ! we need to keep a list of these sites so we can add them
  ! back in later when we print out the list.
    jInactive = jInactive + 1
    inactives(jInactive,:) = (/iD,labelFull(1,iD)/)
  endif
enddo

END SUBROUTINE make_inactive_table

SUBROUTINE getSpaceGroup_activeSitesOnly(pLV,nDfull,dFull,inactives,latDim,rot,shift,eps)
real(dp), intent(in)  :: pLV(3,3)
integer, intent(in)   :: nDFull
real(dp), allocatable :: dFull(:,:)
integer, allocatable  :: inactives(:,:)
integer, intent(in)  :: latDim ! Lat Dimensionality, i.e., surface or bulk
real(dp), allocatable:: rot(:,:,:), shift(:,:)
real(dp), intent(in):: eps

integer, allocatable :: aTyp(:)
integer :: status

allocate(aTyp(nDfull),STAT=status)
if(status/=0)stop "Allocation failed in get_dvector_permutations: aTyp"

aTyp = 1
aTyp(inactives(:,1)) = 2 ! inactive sites should never map to active, so give them a distinct label
call get_spaceGroup(pLV,aTyp,dFull,rot,shift,.false.,eps)
if(latDim==2) call rm_3D_operations(pLV,rot,shift,eps)

END SUBROUTINE getSpaceGroup_activeSitesOnly


END MODULE
