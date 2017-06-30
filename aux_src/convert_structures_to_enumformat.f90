!!<summary>THIS IS ONLY A TEST PROGRAM
!! for converting a structures.in file to structures.enum.out or sth like that.
!! This file is quick and dirty. Copied lots of code from uncle since the steps that the
!! actual program takes will eventually be moved to figure_enumeration</summary>

MODULE ce_types
implicit none

integer, parameter:: dp=selected_real_kind(15,307)
integer, parameter:: sp=selected_real_kind(6,37)
integer, parameter:: si=selected_int_kind(1) ! very short integer -10..10 range
integer, parameter:: li=selected_int_kind(18) ! Big integer -10^18..10^18 range   

type label_pt ! A structure for defining points that become vertices in a figure
   integer  :: label  ! basis point of the parent lattice for this R
   real(dp) :: pos(3) ! Cartesian coordinates of the point
   real(dp) :: len    ! Distance of the point from the origin
endtype label_pt

type label_pt_list
   integer                      :: thisLabel   
   type(label_pt)     , pointer :: thisRList(:)
   integer                      :: thisRListSize
   type(label_pt_list), pointer :: prev, next
end type label_pt_list

type figure ! Stores a CE figure and associated information
   integer  :: nV                   ! number of vertices
   real(dp) :: avgR                 ! average distance of each point from the center of gravity

   real(dp), pointer :: vertex(:,:)  ! { coordinate# , vertex# } 3xn_(Vertex) array of points
   integer, pointer  :: gPt(:,:)     ! { gcoordinate#, vertex# }
                                     ! (4xn_V) the vertex points in
                                     ! gTab coordinates
   integer, pointer  :: label(:)     ! { vertex# } basis point of the
                                     ! parent lattice for each vertex
   integer, pointer  :: s(:)         ! { vertex# } chebychev index for
                                     ! a vertex ("s-vector", when
                                     ! taken for all vertices) See the
                                     ! first couple of pages of PRB
                                     ! 49, 8627 (1993)
   type(figure), pointer :: reference(:)
endtype figure

type figRep  ! Representative figure (with all its sym-equivalent brothers)
   type(figure), pointer :: Rep(:) ! Store all the brothers together
   integer(li) :: count ! Number of representatives (just size of Rep)
                        ! BEWARE: for the pizza figures it is NOT the
                        ! size of represenations but the size of
                        ! equivalent figures (take a look at
                        ! generate_multilattice_pizza
   integer(li) :: norm !  Normalization (Ncells x Nrep figs)
   type(figRep), pointer :: reference(:)
endtype figRep

type t_pizza
   integer               :: nCl
   type(figure), pointer :: cluster(:)
end type t_pizza

type mbeci ! List of ECI's for pairs or higher-order figures
   real(dp), pointer :: J(:)
endtype mbeci

type correlation
   real(dp), pointer :: s(:)  ! each correlation now consists of up 
                              ! to (CErank-1)^(n_Vertex) values, depending on
                              ! which s-Vector of the s-Matrix is applied
                              ! calling this "s" might be a bad idea, but for now ...
endtype correlation

! This structure store various information about each input structure,
! including the results of the CE prediction and ground state search
type crystal
   real(dp) :: LV(3,3)           ! Lattice vectors
   real(dp), pointer :: pos(:,:) ! position of the basis atoms
   integer, pointer :: spin(:)   ! occupation variables (= -1, +1, etc...)
   integer, pointer :: pBVlab(:) ! parent basis vector label
   integer          :: nAtoms    ! total number of atoms
   integer, pointer :: nTyp(:,:) ! The total number of each (spin)
                                 ! type of atom (index: 1 ~ 1st atomic
                                 ! species, 2, 3...)

   logical          :: atomsMissing
   integer          :: missingTyp(2), missingN

   integer, pointer :: xnTyp(:,:)! The total number of each (spin)
                                 ! type of atom. Only the x-relevant
                                 ! atoms are counted
   integer, pointer :: ATyp(:)   ! The (spin) type of each atom
                                 ! (=rank) in the structure = 1, 2,
                                 ! 3...  (1-to-1 mapping to the
                                 ! spin(:))
   character(80) :: title        ! Name of the structure
   real(dp) :: pLV(3,3) ! parent lattice vectors
   real(dp), pointer :: intPt (:,:) ! Interior points of the parent lattice
   real(dp) :: a0                ! Lattice constant (not really necessary but...)
   integer  :: weight            ! weight of the structure during the fit
   real(dp) :: energy            ! Energy of the structure (DFT)
   real(dp) :: ce_energy         ! Energy predicted by CE
   real(dp) :: ce_refEnergy      ! Total Reference energy as provided
                                 ! in all reference structures
   real(dp) :: fitEnergy         ! the energy that will be used for fitting
   real(dp) :: formationenergy   ! Enthalpy of formation of the structure
   real(dp) :: ce_formationenergy! Enthalpy of formation of the structure as predicted by CE
   real(dp) :: quasibinform      ! Quasibinary formation enthalpy, only used in ternaries
   real(dp) :: gstenergy         ! Energy of the groundstate at this
                                 ! concentration (or the line?)
   integer ::  gstleft, gstright ! index in dftgstlist(:) of groundstates to left and right
   real(dp), pointer :: PI(:)
   real(dp) :: eps ! epsilon for finite precision checking
   real(dp) :: conc ! binary concentration, in ternary case: quasibinary conc.
   real(dp) :: con(3) ! ternary concentrations
   real(dp), pointer :: x(:,:)  ! { atomic type, CE# }: concentration
                                ! of each component type (will
                                ! eventually supercede conc and con(3)
                                ! Note: CE# is for coupled CEs. For a
                                ! "normal" CE, the second index is
                                ! always ==1
   real(dp), pointer :: delta(:,:) ! Garbulsky Ceder Constraints: {
                                   ! constraint number, observable }
   logical :: isHelpStructure
   type(crystal), pointer :: reference(:)
   logical                :: automaticReference

   integer, pointer :: HNFlist(:,:,:) ! { entry, entry, HNF# }
   integer, pointer :: pLabel(:,:)    ! { labeling#, entry }

endtype crystal

type surface_variables
  logical  :: isSurf             ! switch: is this a surface CE?

  real(dp) :: unitNormal(3)      ! surface normal vector (user input, then normalised)
  real(dp) :: distT, distB       ! distance of "top" and "bottom"
                                 ! surface from origin (calculated)

!   character(1), pointer :: BasType(:)      ! type of the basis point (bulk "B"/surface "S"/adsorbate "A")
!   integer, pointer      :: BasLayer(:)     ! in which layer is this basis atom?
!   integer, pointer      :: BasLayerAtom(:) ! which atom is it in this layer?
!   integer, pointer      :: BasId(:)        ! lookup table for the site identities of the slab, e.g. atom 1 ~ atom 10 due
!                                            ! to the symmetric slab (and those site always have the same occupancy)
!                                            ! (this does imply that those sites have the same onsite correlations,
!                                            ! see CE_variables%BasEq. However, it does not imply that the sites
!                                            ! with the same onsite correlations, see CE_variables%BasEq, are
!                                            ! identical in the surface slab sense)

!  integer, pointer :: BulkList(:)  ! { bulk atom index, i.e. (1) is the first bulk atom }: which basis atoms are bulk like?
!  integer, pointer :: BulkSpin(:)  ! { dset#, i.e. (1) is the first dset member }: spin variables of the bulk atoms
!  character(1), pointer :: BulkEnumRank(:) ! { dset# }: rank corresponding to spin as a character, enum counting
!                                           ! (i.e. starting at 0 for the first species)

  integer  :: LVpar2normal        ! which surface LV is parallel to the surface normal
!  logical  :: isSymmetric         ! is the slab symmetric?
!  real(dp) :: distSymPlane        ! distance of symmetry plane
                                 
!  real(dp) :: ShiftVec(3)         ! surface shift vector
                                 
!  integer  :: nSlabLayers         ! number of slab layers
!  integer  :: nBulkLayers         ! number of bulk layers
!  integer  :: nSurfLayers         ! number of surface layers

endtype surface_variables

type adsorbate_variables
  logical :: isAds
  logical, pointer  :: isAdsD(:)   ! { dset# }: which of the dset members are adsorbate places?
  
  integer           :: freeSiteRank     ! the rank of the free adsorbate site
  real(dp), pointer :: singleEnergy(:)  ! { (coupling) rank# }:
                                        ! energies of the single
                                        ! adsorbate atoms
!  real(dp), pointer :: cleanEnergy(:)   ! { rank# }: energies of the clean 100% surfaces

end type adsorbate_variables

type lat2ref
  real(dp)           :: Mlv(3,3)           ! transformation matrix for the lattice vectors
  real(dp)           :: Mds(3,3), Ods(3)   ! transformation matrix and
                                           ! offset for the dset
                                           ! members
  integer            :: nDel
  integer, pointer   :: dsetDel(:)
end type lat2ref

type t_individual ! stores fitting parameters, cv score, etc, for one individual
   logical, pointer :: genome(:)  ! T = fitting parameter is activated, F = not activated
   integer, pointer :: onGenes(:) ! only the (indices of the) genes that are switched ON
   integer          :: nOnGenes   ! number of entries in the onGenes list
   real(dp), pointer:: J(:,:)     ! { cluster#, subset# } effective cluster interactions
   logical, pointer :: subsetFitSuccessful(:) ! { subset# }, not really used
   real(dp) :: cvs
   integer  :: age
end type t_individual

type t_fitResults
  integer                     :: bestIndividualIndex 
  type(t_individual)          :: bestIndividual

  real(dp), pointer           :: Jav(:)     ! { cluster# } effective
                                            ! cluster interactions
                                            ! (averaged)
  real(dp), pointer           :: Jff(:)     ! { cluster# } effective
                                            ! cluster interactions
                                            ! (full fit, i.e. all
                                            ! structures, no subsets)

  real(dp), pointer           :: fittingError(:)     ! { structure# }
  real(dp)                    :: fittingErrorMax, fittingErrorAv, fittingErrorRMS
  real(dp), pointer           :: predictionError(:,:)  ! { structure#, subset# }
  real(dp)                    :: CVS
end type t_fitResults

type CE_variables             ! This data type will contain all lat.in data
                              ! and will simplify shoving data back and
                              ! forth information between routines
                              ! a lot more general data should be put 
                              ! in there later for convenience and clarity 
   character(30) :: projectname ! Title of CE-project
   character(10) :: title     ! Title of lattice
   character(1) :: CEtype     ! Bulk or surface?
   logical :: scalar          ! Do you need an energy expansion
                              ! (scalar) only, or something more
                              ! complicated (multidimensional)?
   integer :: nObservables    ! number of configuration dependant
                              ! observables, that is, _including_
                              ! energy
   logical      :: isBulk 
   character(1) :: CoordSys   ! parent basis vectors / surface normals
                              ! in Cart. or Direct coords?
   integer :: CErank          ! Binary, ternary, ... ?
   real(dp) :: pLV(3,3)       ! parent lattice vectors

   type(surface_variables)   :: surface   ! surface variables
   type(adsorbate_variables) :: adsorbate ! adsorbate variables

   integer           :: nBas          ! number of parent basis vectors
   integer           :: nBasIneq      ! number of inequivalent basis vectors / onsitecorrs
   integer,  pointer :: digit(:), label(:,:), equivalencies(:)
   integer,  pointer :: BasEq(:)      ! lookup table showing which
                                      ! basis corresponds to which
                                      ! onsite
   real(dp), pointer :: pBas(:,:)     ! parent basis vectors (3xnBas)
   logical,  pointer :: xRelevant(:)  ! specifies if a dset member is relevant
                                      ! for calculating the concentration

!   integer           :: nBasEnum      ! number of parent basis vectors that will be enumerated
!   logical,  pointer :: BasEnum(:)    ! enumerate this basis point?
!   real(dp), pointer :: pBasEnum(:,:) ! the parent basis vectors for which BasEnum=.true.

   real(dp) :: eps            ! accuracy for handling of coordinates
   real(dp) :: xeps= 0.00000001_dp   ! accuracy for handling of concentrations
   real(dp) :: cutoff(5)      ! cutoffs for figure definition

   logical  :: GSS_override      ! go on even if struct_enum.out does
                                 ! not seem to be consistent with CE
                                 ! data
   logical  :: GSS_concOnly      ! perform full GSS search or only
                                 ! such meeting the required
                                 ! concentrations?
   integer  :: GSS_concN         ! element number for the next variable:
   real(dp) :: GSS_concRange(2)  ! min/max concentration of element GSS_concN in GSS

   logical :: GSS_full        ! perform full GSS search or only "true" multi-naries?
   logical :: GSS_B2          ! perform GSS within a 54-atom B2-like cell?
   integer :: GSS_B2mode      ! evaluated Antisites or Vacancies in B2?
   integer :: GSS_B2defects   ! how many defects allowed in B2 GSS mode?
   integer :: GSSspan(2)      ! Starting and ending cell sizes for the ground state search (GSS)
   logical :: GSS_MPIcat      ! after the GSS: do you wan to
                              ! concatenate the different MPI
                              ! gss.out.* files to a single gss.out?
   logical :: GSS_MPIdelete   ! after the concatenation (see above):
                              ! delete original gss.out.* files?
   
   integer :: totalstructures ! how many structures are defined in struct.in
   integer :: leftoutStr      ! how many structures at the end of structures.in
                              ! are not used for fitting?
   integer, pointer :: pureList(:) ! the pure element structures

   integer, pointer :: GSlist(:)   ! the structures that are ground states
   logical :: quasibin        ! is the ternary CE quasibinary?
   logical :: cstrain, locstrain
   logical :: concdep         ! concentration dependent ECI?
   integer :: polynome        ! polynomial order of ECI x-dependency -1 (!)
   integer, pointer :: dftgstlist(:) ! indices of the groundstate structures

   real(dp), pointer :: delta(:,:,:)             ! { constraint# ,
                                                 ! update# ,
                                                 ! observable# }:
                                                 ! Garbulsky-Ceder
                                                 ! deltas.
                                                 ! constraint# : 1,2,3
                                                 ! for the 3 Deltas
                                                 ! update# : we can
                                                 ! supply a list of
                                                 ! Deltas in fitpar.in
                                                 ! observable# : 1st
                                                 ! is energy
   integer, pointer  :: gen4DeltaSwitch(:)       ! switch to delta*new after x Generations

   real(dp) :: pureDelta1     ! GC-constraint 1 for the pure elements (i.e. the total deviation)
   real(dp) :: gstweight      ! weight of groundstates when performing fit
   real(dp) :: ktfit          ! fit parameter
   real(dp) :: cone, ctwo, lambda1, lambda2      ! fit parameters
   
   real(dp) :: qp_eps         ! epsilon for quadprog routine
   real(dp) :: inf            ! infinity for quadprog routine
   integer  :: maxcounts      ! maxcounts for quadprog routine


   type(figure), pointer       :: fig(:) => null()
   type(figRep), pointer       :: eqvClusters(:) => null()
   type(figRep), pointer       :: eqvClustersList(:,:) => null()
   real(dp), pointer           :: J(:), Jff(:)

   type(CE_variables), pointer :: reference(:)
   integer                     :: nreferences
   character(80)               :: lat_file, J_file, cluster_file, str_file
   type(lat2ref), pointer      :: l2r(:)

   integer                     :: ncouplings   ! how many CEs will be coupled?
   integer, pointer            :: CEbas(:)     ! { dset# }: to what CE do the dvectors belong to?

endtype CE_variables

type t_subset
  integer            :: nFit, nPred
  integer,   pointer :: fitList(:)       ! structures in subset (fit)
  integer,   pointer :: predList(:)      ! structures that are predicted
  integer            :: nGS              ! number of ground states
  integer,   pointer :: GSList(:)        ! ground states
  integer,   pointer :: pureList(:)      ! { element } the indices of the pure elements
  real(dp),  pointer :: GSLDist(:)       ! distance to GS convex hull
  integer,   pointer :: GSofPlane(:,:)   ! { GS idx, structure } the GS that span the hyperplane
                                         ! above which the structure is floating. 
  real(dp),  pointer :: GSMixture(:,:)   ! { mixture idx, structure }
  integer,   pointer :: StrLowE(:)       ! { structure } the idx of
                                         ! the str that is lowest in
                                         ! energy below the current
                                         ! struc
  real(dp),  pointer :: StrLowEDist(:)   ! { structure } distance to
                                         ! structure lowest in energy
                                         ! (fitenergy)
  real(dp),  pointer :: StrDeltaDist(:)  ! { structure } distance to
                                         ! structure lowest in energy
                                         ! (total energy), used for GC
                                         ! deltas
  real(dp),  pointer :: Delta(:,:,:)     ! { delta idx, structure, observable }
  
end type t_subset

type t_GA_continue
  logical       :: yes
  character(80) :: filename
  logical       :: best
  logical       :: dead
  logical       :: all
  logical       :: indv
  integer       :: iIndv
end type t_GA_continue

type GA_variables             ! This data type will contain all GApar.in data
                              ! Not sure if it makes sense NOT to put this into
                              ! CE_variables, but I like it better right now.

   logical :: simplefit       ! perform simple fit or apply GA?

   type(t_GA_continue) :: continue        ! continue a GA run?

   integer :: nClusters, nAlwaysOn        ! total number of clusters,
                                          ! number of clusters that
                                          ! are always on (i.e. the
                                          ! Constant and the Onsites,
                                          ! e.g. monolattice 1+1=2)
   integer :: nMinClusters, nMaxClusters  ! minimum and maximum number of clusters 
                              ! (i.e. figures + s-vectors, i.e. fitting variables!)
   integer :: nTotClusters    ! total number of clusters available
   integer :: popsize         ! size of GA-population
   integer :: maximumAgeForIndividual  
   integer :: nChildren       ! nr. of children generated in each generation
   integer :: nReplacements   ! nr. of children to replace parents
   real(dp):: diversityLimit
   integer :: maxgenerations  ! nr. of generations in each GA run
   integer :: nSubsets        ! nr. of subsets for CVS
   integer :: nPredictionOfEachStructure  ! how often should each structure be predicted?
   integer :: rndseed         ! seed for random number generator
   integer :: fixedseed=99    ! seed from input instead from system clock
                              ! 99=clock, all others: fixed at this value
                              ! Important: this is for OUR RNG ONLY, for MC
                              ! we use the fortran intrinsic random_number
   logical :: includeend = .false.  ! include endpoints of conc-range
                                    ! in predictions? (deprecated,
                                    ! must be false)

   real(dp) :: mutationProbability

   type(t_subset), pointer :: subset(:)   ! prediction/fit subsets

endtype GA_variables


type CE_results             ! This data type will contain all final and temporary
                            ! results of the CE
   real(dp) :: CVS                     ! resulting CVS
   real(dp) :: EVS                     ! resulting EVS
   real(dp), pointer :: ecival(:)  ! the ECI resulting from fit
   real(dp) :: maxerr, averror, rms    ! some errors
   real(dp) :: maxerrpred, averrorpred, rmspred ! more errors for predictions
   real(dp), pointer :: jconst(:)        ! 
   real(dp), pointer :: jcorr(:,:)         ! Just nicer access to final eci
   real(dp), pointer :: mbeci(:,:)       !
   real(dp), pointer :: jconst_first(:)  !
   real(dp), pointer :: jcorr_first(:,:)   !
   real(dp), pointer :: mbeci_first(:,:) !
   integer, pointer :: bestfitset(:)     ! array containing the oprimized figureset
endtype CE_results

type MC_variables ! variables for Monte Carlo
   ! Some general constants first
   integer :: maxBoxSizeForCalcTotalE             = 50
   integer :: time2printSpeedInfo                 = 30
   integer :: timeAvgStep2printProgressBar        = 600
   integer :: SIsteps                             = 100000

   ! MC variables 
   integer :: boxSize(3)                      ! Size of the
                                              ! simulations cell
                                              ! (times unit cell in 3
                                              ! directions)
   character(30) :: cboxSize                  ! character equivalent of boxSize (for printing)
   real(dp), pointer :: x(:)                  ! Concentration of each
                                              ! atom type in the
                                              ! simulation E! NOTE:
                                              ! have to generalise the
                                              ! concentration and mu
                                              ! for coupled CEs
   real(dp), pointer :: mu(:)                 ! chemical potential for each atom type
   real(dp), pointer :: Dmu(:,:)              ! Difference in chemical potential (lookup-table)
   integer(si), pointer :: possibleNewSpins(:,:)  ! { index, current
                                                  ! spin }: list (1st
                                                  ! index) of possible
                                                  ! new spin
                                                  ! occupations if
                                                  ! current spin is
                                                  ! given (2nd index)
                                                  ! e.g.
                                                  ! possibleNewSpins(:,-1)
                                                  ! = (/ 0,1 /) in a
                                                  ! ternary case
   real(dp) :: eps              ! epsilon for finite precision equality checking

   integer :: maxgFigLen        ! maximal length of an interaction,
                                ! measured in gspace (used for
                                ! parallel MC dependency markings)

   integer              :: nTemperatureSchedule
   integer(si), pointer :: TemperatureMode(:) ! 0 = ramp for T,
                                              ! i.e. use temperature
                                              ! decrease factor 1 =
                                              ! Delta for T, i.e. go
                                              ! down linearly
   real(dp), pointer    :: changeT(:)
   real(dp), pointer    :: startT(:), endT(:) ! Start and end temperature for MC-Run

   real(dp) :: E_conv           ! Energy criterion for E(T) convergence
   integer(si) :: nDepLists     ! number of dependency lists
   integer     :: depBoxLength(3)  ! coarse graining size of the dependency marking

   logical     :: onlyTestSpeed ! do we only want to test the speed (and then abort)?
   real(dp)    :: SIstepsFactor

   integer(li)  :: maxnSteps ! maximal number of MC steps
   real(dp), pointer :: norm_vec(:,:) ! (3xn) normal vectors for MCcell evaluation
   real(dp), pointer :: basis_vec(:,:) ! (3xn) basis vectors for MCcell evaluation
   integer  :: planes ! how many planes should be evaluated from MCcell?
   integer  :: rank ! Number of atom types (binary, ternary, etc.)

   integer  :: d(4) ! Extent of the g table in each dimension
   integer  :: D33(3,3)  ! full D matrix (not only the diagonal elements as in d(1:3))

   integer(li)  :: nAt  ! Number of atoms in the sim. cell, prod(d)
   integer, pointer, dimension(:)  :: rk, sp ! Rank space and spin
                                             ! space mapping
                                             ! (-k/2..k/2, 1..k)
   integer :: nBasIneq
   integer,  pointer :: BasEq(:)      ! lookup table showing which
                                      ! basis corresponds to which
                                      ! onsite
   integer(li), pointer :: onsiteNorm(:) ! Number of occurrences of
                                         ! each type of atom (for
                                         ! onsite normalization)

   integer(li) :: avsteps   ! over how many MC-steps will be averaged?
   integer :: nTsteps_write_MCcell ! after how many temperature steps
                                   ! shall a MCcell file be written?
   character(30) :: fname

   logical       :: continue                 ! Continue a previous MC-run?
   character(30) :: continue_FileName        ! MCcell-file to read in 
   integer       :: continue_TemperatureStep ! stepnumber at which calculation will be continued
   real(dp)      :: continue_Temperature     ! Temperature at which continued MC-cell is read in
   real(dp)      :: continue_Energy
   logical       :: continue_EnergyConverged
   integer       :: continue_iTemperatureSchedule
   integer       :: continue_scale(3)        ! Scale factor of the MC
                                             ! continue run (e.g. 2 2
                                             ! 2 => MC cell will be
                                             ! twice the sice of
                                             ! continue cell in each
                                             ! direction

   logical :: GC        ! GrandCanonical Ensemble or Canonical?
endtype MC_variables

type structure_prediction_t
character(1)              :: coordsys       ! coordinate system used for real space output
real(dp), pointer         :: xi(:,:), xf(:,:)  ! concentration ranges: xi = start, xf = stop
logical                   :: automatic      ! automatic search?
integer                   :: nSPX           ! number of structures per
                                            ! concentration (for
                                            ! automatic search)
integer                   :: nS             ! number of structures (for manual search)
integer, pointer          :: struclist(:)   ! the list of enum-numbers of the structures
real(dp), pointer         :: energylist(:)  ! the energy (DHf) of the structures
type(crystal), pointer    :: strucs(:)      ! the structures itself
                                            ! (will be built later
                                            ! than struclist and
                                            ! energylist, using the
                                            ! information in those
                                            ! lists.
integer, pointer          :: GSList(:)      ! the list of groundstates
end type structure_prediction_t


! TK @ GH: can we remove this again???
!
!
!CONTAINS
! Don't call this routine, that is, don't deallocate things that 
! have been copied *from*. The copied things just contain pointers to the original elements. So if
!  ! you deallocate the elements, those copied structures will become corrupted.
!subroutine deallocate_figure_array(fig)
!type(figure), pointer :: fig(:)
!integer i
!do i = 1, size(fig)
!   print *,"deallocating fig #",i
!   call deallocate_figure_elements(fig(i))
!enddo
!deallocate(fig)
!end subroutine deallocate_figure_array
!subroutine deallocate_figure_elements(fig)
!type(figure) fig
!if(associated(fig%vertex))deallocate(fig%vertex)
!if(associated(fig%label))deallocate(fig%label)
!if(associated(fig%gPt))deallocate(fig%gPt)
!if(associated(fig%S))deallocate(fig%S)
!end subroutine deallocate_figure_elements


END MODULE

!!<summary>! test driver file for converting a given structures.in
!! file to enumformat with the help of lat.in the main concept is
!! taken from find_structure_in_list. </summary>
program convert_structures_to_enumformat
use enumeration_utilities
use ce_types
use vector_matrix_utilities
use symmetry
implicit none

type(CE_variables) :: CE
type(crystal), pointer :: str(:)
logical :: match
character(80) POSCAR, LATIN, dummy, title

integer :: SNF(3,3), L(3,3), LatDim ! L is the left transform for SNF

integer iStr, nStr, jStr

call getarg(1,dummy)
read(dummy,'(a80)') LATIN
call getarg(2,dummy)
read(dummy,'(a80)') POSCAR

print *, "LATIN =", LATIN
print *, "POSCAR=", POSCAR

CE%lat_file = LATIN
call read_lattdef(CE,0)
if (CE%surface%isSurf) then; LatDim=2; else; LatDim=3; endif;

call read_input_structures(POSCAR,str,CE)
call check_input_structures(str, CE)
nStr = size(str)

do iStr=1,nStr
  print *, "get HNF of structure #",istr
  call get_HNF_of_derivative_structure(str(iStr)%title,str(iStr)%LV,str(iStr)%pos,str(iStr)%aTyp-1,CE%pLV,CE%pBas,str(iStr)%HNFlist,SNF,L,CE%eps)
  call get_gspace_representation(CE%pLV,CE%pBas,str(iStr)%LV,str(iStr)%pos,str(iStr)%aTyp-1,str(iStr)%HNFlist,LatDim,str(iStr)%pLabel,CE%eps)
enddo

do iStr=1,nStr
do jStr=iStr+1,nStr
  write(*,'(A,I4,I4)',advance='no'), "comparing strucs", istr, jstr
  call compare_two_gstructures(LatDim,CE%pLV,CE%pBas,CE%eps, &
       str(iStr)%HNFlist(:,:,1), str(iStr)%pLabel(1,:), &
       str(jStr)%HNFlist(:,:,:), str(jStr)%pLabel(:,:), &
       match) 
  write(*,'(5x,L)') match
enddo
enddo

contains

  !!<summary>subroutine was taken from the code of Ralf Drautz co_ca
  !!cares about comments and blanks and avoids reading them, comment
  !!lines start with a #</summary>
  subroutine co_ca(unit,error)
    implicit none
    character(50) phrase !letter: contains the first letter of every line
    logical   com !true if comment line is found
    integer  unit, i, ios !unit specifies the unit to read from
    logical error
    
    com = .true.; error = .false.
    do while ( com .eqv. .true.)
       read(unit,50,iostat=ios) phrase
       if (ios/=0) exit
       !              blank line ?
       if (phrase .ne. ' ') then
          i = index(phrase, '#')
          !                 # not found ?
          if (i .ne. 0) then
             !                    # first letter in line ?
             if (i .ne. 1) then
                phrase = phrase(1:i-1)
                if (phrase .ne. ' ') com= .false.
             endif
          else
             com = .false.
          endif
       endif
    end do
50  format(50a)
    backspace unit
    if (ios/=0) error = .true.
  end subroutine co_ca
  
  !!<summary>Makes a string uppercase.</summary>
  subroutine ucase(string)
    implicit none
    character(len=*), intent(INOUT) :: string
    integer :: i
    integer :: length
    
    length=len(string)
    do i=1,length
       if (lge(string(i:i),'a').and.(lle(string(i:i),'z'))) then
          string(i:i)=achar(iachar(string(i:i))-32)
       endif
    enddo
  end subroutine ucase
  
  !!<summary></summary>
  subroutine parse_line_for_numbers(line,separator,maxEntries,readEntries,values)
    character(1000), intent(inout):: line
    character(1), intent(in)      :: separator
    integer, intent(in)           :: maxEntries
    integer, intent(out)          :: readEntries
    real(dp), pointer             :: values(:) ! intent(OUT)                               

    integer iE, ierr

    if (associated(values)) deallocate(values)
    allocate(values(maxEntries))

    readEntries=0
    iE = 0

    line = adjustl(line) ! Remove preceding blanks                                         
    if (trim(line)=="") then ! make sure we didn't get a blank line
       readEntries = 0
       return
    endif
    
    line(1000:1000) = separator

    do while (iE<maxEntries)
       iE = iE+1
       read(line,*,iostat=ierr) values(iE)
       if (ierr /=0) then
          readEntries=iE-1
          return 
          !stop "ERROR: line parsing for number failed. Maybe format mismatch?"
       endif
       line = adjustl(line(index(line,separator)+1:))  ! remove what we just read in
       if (line=="") exit
    enddo
    
    readEntries=iE
  end subroutine parse_line_for_numbers

  !!<summary>This routine takes the CE rank (given by k) and maps it
  !!onto the spin variables, -k/2..k/2, skipping 0 where necessary
  !!(for even rank CEs), and maps it onto rank space, 1..k, as
  !!well</summary>
  SUBROUTINE get_SpinRank_mapping(k, spinSpace, rankSpace)
    integer, intent(in) :: k
    integer, pointer :: spinSpace(:), rankSpace(:)

    integer ik, i ! Counters
    
    allocate(spinSpace(k))
    allocate(rankSpace(-(k-mod(k,2))/2:(k-mod(k,2))/2))  !Handle k=odd and k=even cases together
    do ik = 0, k-1
       i =  ik-k/2 ! convert 0..k-1 label to spin variable -k/2..k/2
       if(mod(k,2)==0) then
          if (.not. i<0) i = i + 1 ! skip 0 as a spin if k is even
       endif
       spinSpace(ik+1) = i 
       rankSpace(i) = iK + 1 ! Add one to make list go 1..k
    enddo
  ENDSUBROUTINE get_SpinRank_mapping

SUBROUTINE read_input_structures(file,str, CE)
character(80) :: file
type(crystal), pointer :: str(:)
type(crystal) :: tStr
!type(crystal) :: tStr(1000)
type(CE_variables), intent(inout) :: CE

real(dp) :: energy
integer iStr, iLV, iL, iAt, c, iTyp, ic, i
integer CErank, indx, ioerr
integer nStr, nAt, nTyp, nSpin, nPos
character(1) :: aBasCoord
character(1000) :: line
logical err, peratom, useweights
integer, pointer :: gettype(:), rankSpace(:)

character(1) :: cTransform
logical :: doTransform
real(dp) :: TransformMatrix(3,3), TransformShift(3)

integer :: nSpecialOptions
integer :: iRef, nRef

real(dp), pointer :: numbers(:)


print *,"Reading in the input structures" 

open(13,file=file,status="old",iostat=ioerr)
if (ioerr/=0) stop "ERROR: Couldn't find 'structures.in' file"
open(14,file="debug_readinstructures.out")
nRef = CE%nReferences

!================================================================================
! Special options for structures.in

nSpecialOptions = 3  ! specify the number here! This is for skipping lines afterwards.

! (1)
call co_ca(13,err)
read(13,*) line ! do you want an expansion for energies only or nonscalar quantities?
call ucase(line)
if(line=="SCALAR") then; CE%scalar = .true.; CE%nObservables = 1;
elseif(line=="MULTIDIMENSIONAL") then 
   CE%scalar = .false.
   call co_ca(13,err)
   nSpecialOptions = nSpecialOptions + 1
   read(13,*) CE%nObservables ! read in the number of additional observables, that is _including_ energy
   write(*,*) "attempting multidimensional cluster-expansion... user beware"
   write(*,*) "structures.in says you are trying to clusterexpand on the free energy and ", CE%nObservables-1, " additional observables"
   ! failsafe
   if (CE%nObservables<2) stop "ERROR: nObservables must be >1."
else; stop "ERROR: in the structures.in file, the first entry should be ""scalar"" or ""multidimensional"""
endif


! (2)
call co_ca(13,err)
read(13,*) line ! are the energies total energies? or given "per atom"?
call ucase(line)
if(line=="PERATOM") then; peratom = .true.
elseif(line=="TOTAL") then; peratom = .false.
else; stop "ERROR: in the structures.in file, the first entry should be ""peratom"" or ""total"""
endif

! (3)
call co_ca(13,err)
read(13,*) line
call ucase(line)
if(line=="WEIGHTS") then; useweights = .true.;
elseif(line=="NOWEIGHTS") then; useweights = .false.
else; stop "ERROR: in the structures.in file, the second entry should be ""weights"" or ""noweights"""
endif
!================================================================================ 

! count how many atom *types* there are
call co_ca(13,err)
do iL = 1,5; read(13,*); call co_ca(13,err); enddo  ! Skip the first 5 lines
read(13,*,iostat=ioerr) cTransform ! transformation wanted?
call ucase(cTransform)
if (cTransform=='T') then
  do iL = 1,4; read(13,*);enddo  ! Skip 4 more lines
else
  backspace(13)
endif
cTransform=""

read(13,'(a1000)') line ! Read in the number of atoms of each type
call parse_line_for_numbers(line," ",maxEntries=CE%CErank*CE%ncouplings,readEntries=c,values=numbers)
CErank = c / CE%ncouplings
write(*,'("structures.in contains structures of CE-rank ",i2)') CErank
if (.not.(CErank==CE%CErank)) then 
   print *, "ERROR: CE rank given in lat.in and structures.in does not match!"
   STOP
end if



! Create a table (gettype) of the atom label for any k-nary
! This is just the list of values that the occupation variables can have
call get_SpinRank_mapping(CErank, gettype, rankSpace)
!c=0
!do iTyp = -CErank/2,CErank/2,+1 !GH We should replace this with a call to get_spinRank_mapping 
!   if (mod(CErank,2)==0 .and. iTyp == 0) cycle
!   c=c+1
!   gettype(c) = iTyp
!enddo
!err = .false.  ! GH What did we need this for?


!--------------------------------------------------------------------------------
! Counting structures
!--------------------------------------------------------------------------------
! start anew
rewind(13)
! We need to skip the special options
do iL=1, nSpecialOptions
  call co_ca(13,err)
  read(13,*) line
enddo


write(*,'(A)') "Counting the structures in structures.in..."

iStr = 0
do ! Read in iStr-th structure   
   call co_ca(13,err)
   if (err) exit ! Jump out of the loop if we are at the end of the file
   iStr = iStr + 1
   write(*,'(I5,2x)',advance='no') iStr
   call read_structure(13,tStr,CE,peratom,useweights,ioerr,isReference=.false.)
   if (ioerr/=0) exit

   allocate(tStr%reference(nRef))
   do iRef=1,nRef
     call read_structure(13,tStr%reference(iRef),CE,peratom,useweights,ioerr,isReference=.true.)
   enddo
enddo ! couting
if (ioerr/=0) then 
  iStr = iStr - 1 ! There was junk at the end of the file---don't count
                  ! it as another input structure.
  print *
  write(*,'(A)') "WARNING: structures.in contains junk input at the end of file."
  write(*,'(A)') "WARNING:   last structure will not be included."
  print *
endif
nStr=iStr
write(*,'("read_input_structures found ",i6," input structures.")') nStr !bch increased to i6
CE%totalstructures = nStr
!--------------------------------------------------------------------------------
! End Counting structures
!--------------------------------------------------------------------------------

allocate(str(nStr))

!--------------------------------------------------------------------------------
! Reading structures
!--------------------------------------------------------------------------------
! start anew
rewind(13)
! We need to skip the special options
do iL=1, nSpecialOptions
  call co_ca(13,err)
  read(13,*) line
enddo
do iStr=1,nStr ! Read in iStr-th structure   
   call co_ca(13,err)
   write(*, '("Reading in structure #",i4, ":  ",$)') iStr

   call read_structure(13,str(iStr),CE,peratom,useweights,ioerr,isReference=.false.)
   if (ioerr/=0) exit

   allocate(str(iStr)%reference(nRef))
   do iRef=1,nRef
     call read_structure(13,str(iStr)%reference(iRef),CE,peratom,useweights,ioerr,isReference=.true.)
!     if (str(iStr)%reference(iRef)%automaticReference) call setup_reference_structure(str(iStr),str(iStr)%reference(iRef),CE)
   enddo

enddo ! reading structures
!--------------------------------------------------------------------------------
! End Reading structures
!--------------------------------------------------------------------------------

close(13)

!--------------------------------------------------------------------------------
! For adsorbates: get single energies
if (CE%adsorbate%isAds) then
!  call read_adsorbateData(CE%CErank,CE%adsorbate%singleEnergy,CE%adsorbate%freeSiteRank)
endif

ENDSUBROUTINE read_input_structures



!================================================================================
! Read a single structure (str) from file (fh)
subroutine read_structure(fh,str,CE,peratom,useweights,ioerr,isReference)
integer, intent(in)  :: fh ! filehandle
type(crystal), intent(out) :: str
type(CE_variables), intent(inout) :: CE
logical, intent(in)  :: peratom, useweights
integer, intent(out) :: ioerr
logical, intent(in)  :: isReference
logical              :: automaticReference
integer iStr, iLV, iL, iAt, c, ic, i
integer CErank, nC, indx
integer nStr, nAt, nSpin, nPos
integer, allocatable :: nTyp(:,:)
character(1) :: aBasCoord
character(80) :: line
logical err

character(1) :: cTransform
logical :: doTransform
real(dp) :: TransformMatrix(3,3), TransformShift(3)

CErank = CE%CErank
nC     = CE%ncouplings



if (.not. isReference) then
! This is NO REFERENCE structure
  call co_ca(fh,err)
  read(fh,'(a80)',iostat=ioerr) str%title  !Read in title
  write(*, '(A50)') adjustl(str%title)
  call co_ca(fh,err)
  read(fh,*) str%a0  ! Read in lattice constant
else
! This IS a REFERENCE structure
  call co_ca(fh,err)
  read(fh,*) line
  call ucase(line)
  if (line(1:1) == 'A') then
    automaticReference = .true.
    str%automaticReference = .true.
  else
    automaticReference = .false.
    str%automaticReference = .false.
    if (line(1:1) /= 'R') then
      print *
      print *, "ERROR: Expected a reference structure following the input structure"
      print *, "ERROR: Reference structure section must start with 'R'"
      stop
    endif
  endif
endif
if (.not. isReference .or. (isReference .and. .not. automaticReference)) then
! If the structure is a REAL structure (and no reference structure)
! OR
! if it is a REFERENCE structure W/O AUTOMATIC determination of the reference structure:
! ==> read inputs
  do iLV = 1, 3 ! Read in lattice vectors
    call co_ca(fh,err)
    read(fh,*,iostat=ioerr) str%LV(:,iLV)
    write(14,'(3(f8.4,1x))') str%LV(:,iLV)
  enddo
  
  
  allocate(str%nTyp(CErank,nC))  ! number of atoms of each type
  allocate(str%xnTyp(CErank,nC)) ! number of atoms of each type: only x-relevant atoms are counted (x=concentration)
           str%xnTyp = 0         ! init
  allocate(str%x(CErank,nC))     ! concentration of each type
     
  ! do we need to transform the lattice vectors (rotate, shift origin?)
  call co_ca(fh,err)
  read(fh,*,iostat=ioerr) cTransform ! transformation wanted?
  call ucase(cTransform)
  if (cTransform=='T') then
    doTransform=.true.               ! yes
    do i=1,3
      call co_ca(fh,err)
      read(fh,*,iostat=ioerr) TransformMatrix(:,i)  ! read transformation matrix
    enddo
    call co_ca(fh,err)
    read(fh,*,iostat=ioerr) TransformShift(:)       ! read transformation (cartesian) shift of atomic positions
                                                    ! i.e. re-set the origin
    str%LV = matmul(TransformMatrix,str%LV)   ! transform the lattice vectors
  else                               ! no transformation
    doTransform=.false.
    backspace(fh)
  endif
  
  
  call co_ca(fh,err)
  read(fh,*) str%nTyp(:,:) ! Read in number of each type of atom for each CE
  write(14,'("Atoms: ",20(i3,1x))') str%nTyp

  
  if (any(str%nTyp<0)) then
    str%atomsMissing = .true.
    if (count(str%nTyp<0) >1 ) stop "ERROR: At most ONE missing atom sort allowed"
    if (str%nTyp(CE%CErank,CE%ncouplings) >= 0) stop "ERROR: Up to now, only the last type can by missing"

    str%missingTyp   = minloc(str%nTyp)
    str%missingN     = abs(str%nTyp( str%missingTyp(1), str%missingTyp(2) ))
    str%nTyp         = abs(str%nTyp)
  else
    str%atomsMissing = .false.
    str%missingN     = 0
  endif

  call co_ca(fh,err)
  read(fh,*) aBasCoord ! Read in coordinate system (Cartesian or direct)
  call ucase(aBasCoord)
  
  nAt = sum(str%nTyp) ! sum over all types of atoms of all CEs = total number of atoms
  str%nAtoms = nAt
  write(14,'("Number of atoms: ",i3)') str%nAtoms

  
  ! Determine concentration (legacy code, should be removed someday)
  !________________________________________________________________________________
str%conc = -10000 ! let's see what happens....
! do we still need this? >  ! Binary CE
! do we still need this? >  if (CErank==2) then ! getting binary concentration
! do we still need this? >    str%conc = real(str%nTyp(1))/real(nAt)
! do we still need this? >  !________________________________________________________________________________
! do we still need this? >  ! Ternary CE
! do we still need this? >  else if (CErank == 3) then ! getting ternary concentrations
! do we still need this? >    str%conc = real(str%nTyp(1))/ &                    ! quasibinary concentration for 
! do we still need this? >         & real((str%nTyp(1)+str%nTyp(3))) ! cases, where 2nd atom type is a vacancy
! do we still need this? >    if ((str%nTyp(1)+str%nTyp(3))==0) str%conc = 0.0
! do we still need this? >    do i=1, 3
! do we still need this? >      str%con(i) = real(str%nTyp(i))/real(nAt)
! do we still need this? >    end do
! do we still need this? >  end if
  !________________________________________________________________________________
  ! Remark: the general concentration str%x(:) is set in check_input_structures.


  ! Read atoms
  !________________________________________________________________________________
  allocate(str%pos(3,nAt))
  do iAt = 1, nAt-str%missingN ! Read in each atom
    call co_ca(fh,err)
    read(fh,*,iostat=ioerr) str%pos(:,iAt)
    if (aBasCoord=='D') then !Convert to lattice coordinates
      str%pos(:,iAt) = matmul(str%LV,str%pos(:,iAt))
      if (doTransform) str%pos(:,iAt) = str%pos(:,iAt) + TransformShift ! if wanted: transform the atomic
                                                                        ! positions (reset origin)
    else ! cartesian coordinates
      if (doTransform) str%pos(:,iAt) = matmul(TransformMatrix,str%pos(:,iAt)) + TransformShift
                                                              ! if wanted: transform by rotating and shifting 
                                                              ! the atomic position (note: in the direct coordinates the
                                                              ! rotation was already included in the lattice vectors, here
                                                              ! in cartesian coordinates we have to rotate the atoms
                                                              ! separately.
    endif
    if (ioerr/=0) then
      write(*,'("Missing atoms in structure ",i4," in ''read_input_structures''")') iStr
      stop
    endif
  enddo
  
  !________________________________________________________________________________
  call co_ca(fh,err)
  
  if (.not. isReference) then ! only read the following information if we are not merely
                              ! reading a reference structure (the energy of which will be
                              ! determined by the code through the use of the reference ECI.
    read(fh,*) str%energy
    
    if (.not. peratom) then ! total energies are given so normalize per atom
      str%energy = str%energy/nAt
    endif
       
    if (useweights) then
      call co_ca(fh,err)
      read(fh,*) str%weight
    else
      str%weight=1
    endif
  endif

  ! Determine atom types and spins for the structure:
  ! set failsafe values, that'll break the code. The correct values are set in check_input_structures.
  allocate(str%ATyp(nAt))
  allocate(str%spin(nAt))
  str%atyp = -100
  str%spin = -100

endif ! no automatic reference

str%isHelpStructure = .false.   ! whenever we actually read in a structure or build a reference structure,
                                ! the structure is "real" and not only used for helping some other stuff (e.g. 
                                ! convex hull)

end subroutine read_structure



!***************************************************************************************************
recursive SUBROUTINE read_lattdef(CE,MPI_rank,iReference) ! title,cetype,CErank,LV,aBas,aBasCoord,Ninterior,cutoff,eps)
use utilities, only: ralloc

type(CE_variables), intent(inout) :: CE
integer, intent(in)               :: MPI_rank
integer, intent(in), optional     :: iReference

logical                           :: isReference
logical pr ! print statments?
integer fh ! filehandle

integer ivec, ipt, ipt1, ipt2, iRef, icoord
logical :: err
real(dp) :: cart2latt(3,3)

integer ioerr
integer iBulkAtoms
integer iTopAtom, iBottomAtom

character(1000) :: line
integer, pointer :: SpinSpace(:), RankSpace(:)

if (present(iReference)) then; isReference=.true.; else; isReference=.false.; endif
if (MPI_rank == 0) then; pr=.true.; else; pr=.false.; endif

fh = 10
if (isReference) fh=fh+iReference

open(fh,file=CE%lat_file,status='old',iostat=ioerr)
if (ioerr /= 0) then
  print *, "ERROR: Opening lattice file ", trim(CE%lat_file)
  stop
endif
if (pr) write(*,'(A,A)') "Reading lattice definition from ", adjustl(trim(CE%lat_file))
print *

!-------------- title ---------------------
call co_ca(fh, err)
read(fh,*) CE%title

!-------------- CE type  ------------------
call co_ca(fh, err)
read(fh,*) CE%cetype ! Surface or bulk or adsorbate?
call ucase(CE%cetype)
CE%isBulk          = .false.
CE%surface%isSurf  = .false.
CE%adsorbate%isAds = .false.
CE%ncouplings      = 1
!----- surface
if (CE%cetype=='S') then ! this is a surface CE
  CE%surface%isSurf=.true.
!----- adsorbate
elseif (CE%cetype=='A') then
  CE%CEtype = 'S'            ! will be used for the enumeration, the fact that it is an adsorbate CE (i.e. coupled
                             ! is taken care of by the next variables)
  CE%surface%isSurf =.true.  ! an adsorbate problem is also a surface problem
  CE%adsorbate%isAds=.true.
  CE%ncouplings     =2       ! the CE will be 2 coupled CEs (IMPORTABERT: the rank of each CE must be
                             ! EQUAL at the moment...)
!----- bulk
elseif (CE%cetype=='B') then
  CE%isBulk = .true.
!----- don't know
else
  stop "ERROR: type of CE undetermined. Must be 'B' or 'S' or 'A'"
endif

!-------------- CE rank -------------------
call co_ca(fh, err)
read(fh,*) CE%CErank
call co_ca(fh, err)
read(fh,*) CE%nObservables
!-------------- epsilon -------------------
call co_ca(fh, err)
read(fh,*) CE%eps
!-------------- lattice vectors -----------
call co_ca(fh, err)
call read_pLV(fh,CE%pLV)
!-------------- basis atoms ---------------
call co_ca(fh, err)
read(fh,*) CE%nBas

!----- surface
!if (CE%surface%isSurf) then
!  call get_SpinRank_mapping(CE%CErank, SpinSpace, RankSpace)
!  allocate(CE%surface%BasType(CE%nBas))
!  allocate(CE%surface%BulkSpin(CE%nBas))
!  allocate(CE%surface%BulkEnumRank(CE%nBas))
!  allocate(CE%surface%BulkList(CE%nBas))
!  CE%BasEnum = .false.  ! the correct enumeration will be determined later
!----- bulk
!else
!  CE%BasEnum = .true.   ! will enumerate all dset members
!endif

call co_ca(fh, err)
read(fh,*) CE%CoordSys
call ucase(CE%CoordSys)


!________________________________________________________________________________
! Read in: dset
call read_dset(fh,CE%CErank,CE%nBas,CE%pBas,CE%digit,CE%label,CE%xRelevant,CE%equivalencies)
call transform_dset(CE%CoordSys,CE%pLV,CE%pBas,CE%eps)   ! e.g. direct -> cartesian

! Make sure that all basis atoms are at different locations
call check_for_collisions(CE%pBas)

! Preliminary CEbas setup 
allocate(CE%CEbas(CE%nBas))
CE%CEbas = 1   ! a priori all d-vectors belong to the 1st CE
!---------- end of basis atoms read ----------

!----- surface
if (CE%surface%isSurf) then ! surface CE?

  if (pr) then
    write(*,*)
    write(*,*) "This is a surface CE"
  endif

  call read_surface_normal(fh,CE%pLV,CE%surface%unitNormal)
  call find_topbottom(CE%pLV,CE%pBas,CE%surface%unitNormal,CE%surface%distT,CE%surface%distB)

  if (pr) then
    write(*,'(1x,A,3(F10.3,1x))') "* Surface normal:", CE%surface%unitNormal
    write(*,*) "* Found surfaces at distances from origin:"
    write(*,'(A4,3(A15,A2))')   "|","Bottom Surface", "|", "Top Surface",   "|"
    write(*,'(A4,3(F15.3,A2))') "|",CE%surface%distB, "|", CE%surface%distT,"|"
    print *
  endif

  ! if adsorbate CE: read which dset members are adsorbate sites
  allocate(CE%adsorbate%isAdsd(CE%nBas))
  if (CE%adsorbate%isAds) then 
    call read_adsorbate_sites(fh,CE%nBas,CE%adsorbate%isAdsd)
    where (CE%adsorbate%isAdsd)
      CE%CEbas = 2 ! the adsorbate d-vectors belong to CE 2
                   ! (maybe we could use this variable ONLY, but for the time being I consider it
                   ! more helpful to have the isAdsd too...)
    end where
  else
    CE%adsorbate%isAdsd = .false.
    CE%CEbas            = 1
  endif

endif ! surface CE
!----- bulk
! do nothing

call co_ca(fh, err)
read(fh,*) CE%cutoff

call co_ca(fh,err)
read(fh,*) CE%nreferences

if (isReference .and. CE%nreferences>0) stop "ERROR: References only allowed in the master lat.in"

allocate(CE%reference(CE%nreferences))
allocate(CE%l2r(CE%nreferences))

do iRef=1,CE%nreferences
  if (pr) write(*,'(A,I5)') "Setting up reference ", iRef
  call co_ca(fh,err)
  read(fh,*) CE%reference(iRef)%lat_file
  call co_ca(fh,err)
  read(fh,*) CE%reference(iRef)%cluster_file
  call co_ca(fh,err)
  read(fh,*) CE%reference(iRef)%J_file

! Not yet implemented. As for now I use adsorbate.in for the relevant information.  
!  ! in adsorbate mode we need the structures.in of the surface reference
!  if (CE%adsorbate%isAds .and. CE%reference(iRef)%surface%isSurf) then
!    call co_ca(fh,err)
!    read(fh,*) CE%reference(iRef)%str_file
!  endif

  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Mlv(:,1)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Mlv(:,2)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Mlv(:,3)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Mds(:,1)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Mds(:,2)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Mds(:,3)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%Ods(:)
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%nDel; allocate(CE%l2r(iRef)%dsetDel(CE%l2r(iRef)%nDel))
  call co_ca(fh,err);  read(fh,*) CE%l2r(iRef)%dsetDel(:)
  call read_lattdef(CE%reference(iRef),MPI_rank,iRef)
enddo

!if (CE%surface%isSurf) then
!  if (CE%nReferences == 0) print *, "WARNING: consider using a bulk reference for surface CE"
!  if (CE%nReferences > 1) stop "ERROR: for surface CE only 1 (bulk) reference can be used now"
!  if (any(CE%reference(:)%CEtype /= 'B')) stop "ERROR: for surface CE only bulk references are implemented"
!endif
  
if (CE%adsorbate%isAds .and. .not. any(CE%reference(:)%surface%isSurf)) & ! adsorbate CE but no surface reference
     stop "ERROR: an (bi-binary) adsorbate CE needs a surface reference!"

!call co_ca(fh,err)
!read(fh,*) CE%GSSspan

!REMARK: This should be read in someday, but for now I always want this to be true here
!CE%GSS_full = .true.

close(fh)

END SUBROUTINE read_lattdef

function ProjectionLength(vec,direction)
  real(dp)  :: ProjectionLength
  real(dp)  :: vec(3)
  real(dp)  :: direction(3)
  
  direction = direction / norm(direction)
  ProjectionLength=abs(dot_product(vec,direction))
  
end function ProjectionLength

subroutine read_pLV(fh,pLV)
  integer, intent(in)   :: fh
  real(dp), intent(out) :: pLV(3,3)
  integer i
  logical err
  
  do i=1,3
    call co_ca(fh, err)
    read(fh,*) pLV(:,i)
  enddo
end subroutine read_pLV



subroutine read_dset(fh,k,nD,d,digit,label,xrel,eq)
  integer, intent(in) :: fh          ! filehandle
  integer, intent(in) :: k           ! rank
  integer, intent(in) :: nD          ! number of dset points
  real(dp), pointer   :: d(:,:)      ! intent(out): { coord, dset# } dset points
  integer, pointer    :: label(:,:)  ! intent(OUT): { label#, dset# } labels for the enumeration
  integer, pointer    :: digit(:)    ! intent(OUT): { dset# } number of labels for each point
  logical, pointer    :: xrel(:)     ! intent(OUT): { dset# } relevant for concentration?
  integer, pointer    :: eq(:)       ! intent(OUT): { dset# } equivalency list for enumeration

  integer, pointer :: xrelTmp(:)
  integer :: nE
  real(dp), pointer :: numbers(:)

  character(1000) :: line
  logical err
  integer iD, icoord, i
  
  allocate(d(3,nD),label(k,nD),digit(nD),eq(nD))

!--------------------------------------------------------------------------------
! loop over dset  
  do iD=1,nD
    
    call co_ca(fh, err)
    
    read(fh,'(a1000)') line
    call parse_line_for_numbers(line," ",maxEntries=3,readEntries=nE,values=numbers)
    if (nE /= 3) stop "ERROR: parsing line with d-vector coordinates in lat.in"
    do icoord=1,3
      d(icoord,iD) = numbers(icoord)
    enddo
    call parse_line_for_numbers(line,"/",maxEntries=k,readEntries=nE,values=numbers)
    if (nE==0) then ! only the d-vector coordinates were given
       digit(iD) = k 
       label(:,iD) = (/ (i,i=1,k) /) -1
    else ! occupation was supplied too
       do i = 1, k ! Loop over the number of (possible) labels, exit when there are no more /'s
         label(i,iD) = nint(numbers(i))
   
         ! Sanity check on the input for the label (perhaps not sufficient but catches some errors)
         if (label(i,iD) > k-1 .or. label(i,iD) < 0) then
            write(*,'("Incorrect number for a label, ''",i2,"'', on d-vector #",i2)') label(i,iD), iD
            stop
         endif
         if (i==nE) exit
       enddo
       digit(iD) = i ! Store the number of labels that were specified for each d-vector
    endif
  enddo ! iD
!--------------------------------------------------------------------------------

  ! read which dset members are relevant for concentration
  allocate(xrel(nD)); xrel=.true.
  where (digit==1)
    xrel=.false.
  end where
!  call co_ca(fh,err)
!  read(fh,'(A100)') line
!  allocate(xrel(nD)); xrel=.false.
!  if (line(1:1)=='C') then
!    call co_ca(fh,err)
!    read(fh,'(A100)') line
!    call parse_line_for_numbers(line," ",maxEntries=nD,readEntries=nE,values=numbers)
!    allocate(xrelTmp(nE))
!    xrelTmp = nint(numbers(:nE))
!    do iD=1,nD; if (any(xrelTmp==iD)) xrel(iD)=.true.; enddo
!  else
!    xrel = .true.
!    backspace(fh)
!  endif

  ! read equivalencies (if present, as denoted by "E")
  call co_ca(fh,err)
  read(fh,'(A1000)') line
  if (line(1:1)=='E') then
    call co_ca(fh,err)
    read(fh,*) eq(:)
  else
    eq = (/(i,i=1,nD)/)
    backspace(fh)
  end if

end subroutine read_dset

subroutine transform_dset(coordsys,LV,d,eps)
  character(1), intent(in) :: coordsys
  real(dp), intent(in) :: LV(3,3)
  real(dp), intent(inout) :: d(:,:)
  real(dp), intent(in) :: eps
  
  real(dp) :: LVinv(3,3)
  logical  :: err
  integer ipt
  
  ! convert direct coordinates to cartesian coordinates
  if (CoordSys=='D') d = matmul(LV,d)
 
  ! If the basis atoms are not inside unit cell 0, then shift by lattice vectors until they are
  call matrix_inverse(LV, LVinv, err)
  do ipt = 1, size(d,2)
    call bring_into_cell(d(:,ipt), LVinv, LV, eps)
  enddo
  
end subroutine transform_dset

subroutine check_for_collisions(d)
  real(dp), intent(in) :: d(:,:)
  integer ipt1, ipt2, nD
  
  nD=size(d,2)

  do ipt1=1, nD
    do ipt2=1, nD
      if (ipt1==ipt2) cycle
      if (norm(CE%pBas(:,ipt1)-CE%pBas(:,ipt2))==0) stop "Error: Basis atoms collision!"
    enddo
  enddo
end subroutine check_for_collisions

subroutine read_surface_normal(fh,LV,normal)
  integer, intent(in)   :: fh         ! filehandle
  real(dp), intent(in)  :: LV(3,3)    ! lattice vectors (parent)
  real(dp), intent(out) :: normal(3)  ! normal
  character(1) :: coordsys
  logical :: err
  
  call co_ca(fh,err)
  read(fh,*) coordsys
  call ucase(coordsys)
  
  call co_ca(fh,err)
  read(fh,*) normal
  
  if (coordsys=='D') normal=matmul(LV,normal)

  normal = normal / norm(normal)   ! make unitnormal

end subroutine read_surface_normal

subroutine read_topbottom(fh,d,n,distT,distB)
  integer, intent(in)   :: fh             ! filehandle
  real(dp), intent(in)   :: d(:,:)        ! dset
  real(dp), intent(in)  :: n(3)           ! normal
  real(dp), intent(out) :: distT, distB   ! distances to top and bottom layer, resp.

  logical :: err
  integer :: iTop, iBottom

  call co_ca(fh,err)
  read(fh,*) iTop, iBottom

  distT = ProjectionLength(d(:,iTop)   , n)
  distB = ProjectionLength(d(:,iBottom), n)

  
end subroutine read_topbottom


! Find the top and the bottom surface distance
! by Sascha Maisel
subroutine find_topbottom(pLV,d,n,SdistT,SdistB)
real(dp), intent(in) :: pLV(3,3)  ! parent lattice vectors
real(dp), intent(in) :: d(:,:)    ! d-vectors
real(dp), intent(in) :: n(3)      ! surface normal (unit normal!)
real(dp), intent(out) :: SdistT, SdistB  ! top and bottom distance of the surfaces

integer :: iD,nD
real(dp) :: a(3), b(3)
real(dp) :: M(3,3) ! tranformation matrix to planar basis
real(dp), allocatable :: dslash(:,:) !transformed coordinates
real(dp), allocatable :: tmp(:),trueDistance(:), internalDistance(:)
real(dp) :: tmp1, tmp2
integer, pointer :: list(:)
integer :: gapParameter

nD = size(d,2)

tmp1 = 0._dp
tmp2 = 0._dp
gapParameter = nD

allocate (list(nD))
allocate (dslash(3,nD))
allocate (trueDistance(nD), tmp(nD), internalDistance(nD))

! reconstruction of planar basis:
! Take n to be one of the coordinate directions, and construct a and b accordingly to form a 
! full coordinate system. The n-direction is the first coordinate.
a = (/ -(n(3)**2+n(2)**2), -n(1)*n(2), n(1)*n(3) /)
a = a*1.0_dp / norm(a)

b = (/ 0._dp, n(3), -n(2) /)
b = b*1.0_dp / norm(b)

!inversion of transformation matrix
M(1,1:3) = n
M(2,1:3) = a
M(3,1:3) = b

!perform rotation 
forall(iD=1:nD)
  dslash(:,iD) = matmul(M, d(:,iD))
end forall
trueDistance(:) = dslash(1,:)   ! remember: the first coordinate in the new system is the n direction

if (nD==1) then ! only one d-vector: the distances of the surface(s) is quite clear. We have to issue
                ! return, since the following algorithm relies on at least 2 d-vectors
  SdistT = trueDistance(1)
  SdistB = trueDistance(1)
  return
endif

! sort the points according to their distance
list = (/(iD,iD=1,nD)/)
!call heapsort(list, trueDistance)  ! trueDistance is now sorted...

!calculate intermediate distances (with respect to the n-direction)
forall(iD=1:nD)
  internalDistance(iD) = abs(trueDistance(iD) - trueDistance(mod(iD, nD) + 1))
end forall

!figure out between which atoms the vacuum is located
do iD=1, nD
   if(internalDistance(gapParameter) < internalDistance(iD)) gapParameter = iD
enddo !after this, we can be sure that the vacuum is between gapParameter and gapParameter+1

if(gapParameter == nD) then !this means: vacuum above and under the system
   sDistT = trueDistance(nD)
   sDistB = trueDistance(1)
else !this means: vacuum inside the system
   sDistT = trueDistance(gapParameter + 1)
   sDistB = trueDistance(gapParameter)
endif
end subroutine find_topbottom


subroutine read_adsorbate_sites(fh,nD,ads)
integer, intent(IN) :: fh      ! filehandle
integer, intent(IN) :: nD      ! number of dset points
logical, pointer    :: ads(:)  ! INTENT(OUT), { dset# }: which of the dvectors are adsorbate sites?

logical err
character(1000)    :: line
integer           :: nE
real(dp), pointer :: numbers(:)

if (.not. associated(ads)) allocate(ads(nD))
ads = .false.  ! a priori: no dvectors are adsorbate sites

call co_ca(fh,err)
read(fh,'(A1000)') line ! read line

call parse_line_for_numbers(line," ",maxEntries=nD,readEntries=nE,values=numbers) ! get the numbers of the adsorbate sites
if (nE==0) stop "ERROR: you're trying an adsorbate CE without any adsorbate"

ads(nint(numbers(:nE))) = .true. ! set "ads" of the adsorbate sites to .true.

end subroutine read_adsorbate_sites



logical function isThereAtom(d,x,n,eps)
! this function finds out if there is an atom in the plane (defined by the surface normal) 
! at point x(:)
real(dp), intent(in) :: d(:,:)  ! d-vectors
real(dp), intent(in) :: x(3)    ! actual position
real(dp), intent(in) :: n(3)    ! surface unit normal
real(dp), intent(in) :: eps     ! finite precision parameter

real(dp)  :: distV(3)  ! distance vector
integer   :: iD, nD

isThereAtom=.false.
nD = size(d,2)

do iD=1,nD
  DistV(:)=x(:)-d(:,iD)
  if (ProjectionLength(DistV,n)<eps) then
    isThereAtom=.true.
    return
  endif
enddo
end function isThereAtom










recursive SUBROUTINE check_input_structures(str, CE, thisRef)
use compare_structures
type(crystal), pointer :: str(:), currstr, refstr(:)
type(crystal)          :: mirrorstr
type(CE_variables), intent(in) :: CE
integer, optional :: thisRef
integer :: iRef, nRef
logical :: isReference

real(dp) :: conVec(3), pVol, currD(3), currPos(3), invLV(3,3)
integer sVol
integer iStr, jStr, iAt, jAt, iuq, dlabel, c, iC, nC, iD, nD
integer nStr, nPar, nAt, x, y, z
integer, allocatable :: uqlist(:)
logical checkPoint, equivalent, found, mapped, err, performedReconstruction

integer idx1, idx2, idxoff, ik, iat1, iat2, i
integer, pointer :: gettype(:), rankSpace(:), missingD(:)

character(80) :: header(4), filename

call get_SpinRank_mapping(CE%CErank, gettype, rankSpace)


nC = CE%ncouplings

if (present(thisRef)) then; isReference=.true.; else; isReference=.false.; endif
if (.not. isReference) then
  write(*,'(A)') "Checking input structures"
endif

nStr = size(str)
if (nStr==0) stop "ERROR: check_input_structures: No input structures" 
nPar = size(CE%pBas,2)

do iStr = 1, nStr

   currStr => str(iStr)
   
   if (currStr%isHelpStructure) cycle  ! do not check mere help structures

   if (.not. isReference) then
     write(*,'(A,1x,I5)') "Checking structure #", iStr
   else
     write(*,'(A,1x,I5,1x,A,1x,I5)') "Checking Structure #", iStr, " Reference #", thisRef
   endif


   ! Is each structure a derivative superlattice of the parent lattice?
   if (.not. is_derivative(CE%pLV,currStr%LV,CE%eps)) then
     if (.not. isReference) then
       write(*,'("ERROR: Structure #: ",i4," is not a derivative superlattice.")') iStr
       write(*,'(A,A)') "ERROR:    struc title: ", trim(currStr%title)
       stop
     else
       write(*,'(A,I5,A)') "ERROR: the reference #",iRef," structure is not a derivative superlattice."
       stop
     endif
   endif

   nAt     = size(currStr%spin) ! total number of atoms
   nD      = CE%nBas  ! number of dvectors
   sVol    = nAt / nD ! volume of this structure in terms of unit cells
   

   currStr%xnTyp = 0

   allocate(currStr%pBVlab(nAt))  ; currStr%pBVlab = -1

   ! Do all the interior points of the input structure correspond to points of the parent?
   do iAt = 1, nAt-currStr%missingN

      ! get the dset label of that atom
!      dlabel = get_dset_label( CE%pLV, currStr%pos(:,iAt), CE%pBas(:,:), CE%eps )
      currStr%pBVlab(iAt)=dlabel

      ! this atom of this structure does not lie on a lattice point:
      if(dlabel<0) then
        if (.not. isReference) then
          write(*,'(A,I4,A)') "ERROR: in input structure #", iStr
        else
          write(*,'(A,I4,A)') "ERROR: the reference #",iRef," structure"
        endif
        write(*,'(A,A80)')    "ERROR: name:  ",  adjustl(currStr%title)
        write(*,'(A,I4,A)')   "ERROR:   atom ", iAt, " is not on the parent lattice."
        write(*,'(A,3F5.2)')  "ERROR:   cartesian coordinates of atom: ", currStr%pos(:,iAt)
        stop
      endif

      ! 
      ! if we are here: this atom lies on a lattice point, all is ok (up to now at least...)
      !

   enddo ! loop over all atoms in this structure

   ! now reconstruct missing atoms:
   ! (note: this is a brute force method, not really elegant, but works...)
   performedReconstruction=.false.
   if (currStr%missingN>0) then
      performedReconstruction=.true.
      if (currStr%missingN>50) then
        write(*,'(A)')   "    ... reconstructing atoms: this can take a while for >50 atoms"
      else
        write(*,'(A)')   "    ... reconstructing atoms"
      endif
      
      allocate(missingD( currStr%missingN )) ! allocating an array that is going to held the missing dvectors.
                                             ! If a dvector is missed more than once, it appears more than once!
      missingD = 0
      
      idx1=1
      do iD=1,nD ! loop through the dset
        idx2 = idx1 + (sVol - count(currStr%pBVlab==iD)) -1
                      ! <----------------------------> !
                              ! this expression tells you how often a specific dvector is missing. It should be
                              ! there the same number of times as we have unit cells (sVol = structure volume).
        if (idx2>=idx1) then  ! If it is missing, we ...
          missingD(idx1:idx2) = iD   ! ... add the dvector the correct number of times to the array
          idx1 = idx2+1
        endif
      enddo
      iAt = nAt-currStr%missingN
      do c=1,currStr%missingN  ! loop through all missing atoms. The following is not really elegant, but works.

        if (currStr%missingN > 50) write(*,'(I5,1x)',advance='no') c
   
        shiftloop: do x=0,sVol; do y=0,sVol; do z=0,sVol;
          currD = CE%pBas(:, missingD(c)) ! take the dvector...
          currD = currD + x*CE%pLV(:,1) + y*CE%pLV(:,2) + z*CE%pLV(:,3) ! ... and shift it by parent lattice vectors
          currStr%pos(:,iAt+c)  = currD         ! add the current dvector coordinate to the structure
          currStr%pBVlab(iAt+c) = missingD(c)   ! and set the label correctly
!          call check_for_duplicate_atoms(currStr, iat1, iat2, CE%eps, limit=iAt+c)   ! check if added basis point
                                                ! collides with an already existing point. If yes, the two integers
                                                ! iat1 and iat2 provide the ID numbers of the atoms that collide
          if (iat1==0 .and. iat2==0) then  ! if no collision: both iat1 and iat2 == 0 and ...
            exit shiftloop ! ... we can go on with the next missing atom.
          endif
        enddo; enddo; enddo shiftloop;

        ! failsafe check:
        if (x==sVol .and. y==sVol .and. z==sVol) stop "ERROR: bad programming in the reconstruction loop"
      enddo
   
   endif ! reconstructing
   if (currStr%missingN > 50) write(*,*)

   ! Determine atom types and spins for the structure
   if (.not. isReference) then ! for references this is done during the buildup of the reference structure itself
     do iAt=1,nAt
       currStr%atyp = 1              ! setu: all atoms have the same type...
       currStr%spin = gettype(1)     ! ... and the respective spin
       do iC = 1, CE%ncouplings
         do ik = 2, CE%CErank       ! now: for all other types...
           
           idxoff = sum(currStr%nTyp(:,:iC-1))                    ! all the atoms that were already counted
           idx1   = sum(currStr%nTyp(1:ik-1,iC))+1  + idxoff 
           idx2   = sum(currStr%nTyp(1:ik,iC))      + idxoff
           
           currStr%atyp( idx1 : idx2 ) = ik            ! ... set the type
           currStr%spin( idx1 : idx2 ) = gettype(ik)   ! ... and spin of those atoms
         enddo
       enddo
     enddo
   endif



   do iAt=1,nAt

      ! check if the occupation of every site is ok with the possible occupations in lat.in
      ! CE%label(:,iD) has the possible occupations labels (rank space) for dset member iD.
      ! Thus, the following checks if any one of the possible labels fits to the structure's occupation.
      ! Careful: %label = 0,1,2,3...,k-1 like in enumlib
      !          %aTyp  = 1,2,3,...k
!      dlabel = currStr%pBVlab(iAt)
!      if (.not. any(CE%label(:,dlabel)==currStr%ATyp(iAt)-1)) then
!        if (.not. isReference) then
!          write(*,'(A,I4,A)')   "ERROR: in input structure #", iStr
!        else
!          write(*,'(A,I4,A)')   "ERROR: in reference #",iRef," structure"
!        endif
!        write(*,'(A,A80)')      "ERROR: name:  ", currStr%title
!        write(*,'(A)')          "ERROR:   the structure's occupation does not fit the occupation in lat.in."
!        write(*,'(A,I3,A,I2)')  "ERROR:   structures.in: structure atom        ", iAt,  " has type  ", currStr%ATyp(iAt)-1
!        write(*,'(A,I3,A,10(I2,1x))') &
!                                "ERROR:   lat.in       : equivalent basis atom ", dlabel," has types ", &
!                                                              CE%label(:CE%digit(iAt),dlabel)
!        stop 
!      endif
            
      ! update xnTyp, which stores the total number of atoms of a certain kind, but only
      ! counts x-relevant (i.e. concentration-relevant) structures
!      if (CE%xRelevant(dlabel)) then ! if this is a concentration-relevant dset member
!        currStr%xnTyp(currStr%ATyp(iAt),CE%CEBas(dlabel))=currStr%xnTyp(currStr%ATyp(iAt),CE%CEBas(dlabel)) + 1
!      endif
   enddo

!OLD>   !____________________________________________________________
!OLD>   ! is this a surface CE with symmetric slab?
!OLD>   if (CE%surface%isSurf) then
!OLD>      if (CE%surface%isSymmetric) then
!OLD>
!OLD>        if (.not. structure_is_symmetric(currStr%pos, currStr%spin, CE%surface, CE%eps)) then
!OLD>          write(*,'(A,I4,A)')     "ERROR: in surface input structure #", iStr
!OLD>          write(*,'(A,A80)')      "ERROR: name:  ", currStr%title
!OLD>          write(*,'(A)')          "ERROR:     the structure is not symmetric"
!OLD>          stop
!OLD>        endif
!OLD>
!OLD>      endif ! symmetric slab
!OLD>    endif ! surface CE
!OLD>   !____________________________________________________________


   ! Store parent lattice information for each input structure
   currStr%pLV = CE%pLV 
   allocate(currStr%intPt(3,size(CE%pBas,2)))
   currStr%intPt = CE%pBas
   currStr%eps = CE%eps

   ! update concentration
   !    calculate the concentration from xnTyp, which contains the number of concentration-relevant
   !    atoms in this structure. So, %x contains as many components as there are atom types.
   do iC=1,nC
     currStr%x(:,iC) = currStr%xnTyp(:,iC) * 1.d0 / sum(currStr%xnTyp(:,iC))
   end do
   !    for compatibility reasons: %conc is the quasi-binary concentration.
   currStr%conc = currStr%x(1,1)  

enddo ! loop over all structures, iStr


if (.not. isReference .and. performedReconstruction) then
   filename  = "structures.out"
   header(1) = "# Postprocessed structures.in file"
   header(2) = "#    This file is only written if some of the original atoms were marked"
   header(3) = "#    missing by the use of -n in structures.in, where n is the number of"
   header(4) = "#    missing atoms"
!   call write_InpStrucList(filename,str,header)
endif


!--------------------------------------------------------------------------------
! (in routine: check_input_structures)
! References:
! Now setup the reference structure for the current input structures, check if
! they are ok with the lattice definition (etc.)
if (.not. isReference) then
  nRef=CE%nreferences
!  if (associated(refStr)) deallocate(refStr)
  allocate(refStr(nStr))
  do iRef=1,nRef ! setup all references...
    do iStr=1,nStr ! ... for all input structures
      if (str(iStr)%reference(iRef)%automaticReference) then
        write(*,'(A,1x,I5,1x,A,1x,I5)') "Setting up automatic reference #", iRef, " for structure #", iStr
!        call setup_reference_structure(str(iStr),refStr(iStr),CE,iRef)
        allocate(refStr(iStr)%nTyp(CE%CErank,CE%ncouplings))   ! don't want to do these allocs in the
        allocate(refStr(iStr)%xnTyp(CE%CErank,CE%ncouplings))  ! setup_reference_structure routine, since it is
        allocate(refStr(iStr)%x(CE%CErank,CE%ncouplings))      ! also going to be used during gss
      else
        write(*,'(A,1x,I5,1x,A,1x,I5)') "Setting up user-input reference #", iRef, " for structure #", iStr
        refStr(iStr) = str(iStr)%reference(iRef)  
      endif
    enddo ! end: setup 1 reference (iRef) for all input structures

    ! check the reference structures for reference iRef:
    call check_input_structures(refStr, CE%reference(iRef),iRef)

    do iStr=1,nStr
      str(iStr)%reference(iRef)=refStr(iStr)
    enddo
  enddo ! iRef
  deallocate(refStr)
endif


! Check that each input structure is unique from every other
if (.not. isReference) then
open(100,file='duplicate_structures.out',status='replace')
write(*,'("Checking the uinqueness of the input structures")')
write(*,'(i4," structures in the input list")') nStr
do iStr = 1, nStr-1
   write (*,'("Checking uniqueness of Str. #",I4,"  ",A50)') iStr, adjustl(str(iStr)%title)
   do jStr = iStr+1, nStr
!      write(*,'("Comparing strs. #: ",2(i3,1x))') iStr, jStr

      equivalent=.false.  ! tk: actually, shouldn't need this, but somethings going wrong in 
                          !     compare_arbitrary_structures, where it is set =.false. if
                          !     this optional argument is present. But it fails quite often. 
                          !     Don't ask ME about it...

      call compare_arbitrary_structures(str(iStr)%LV,str(iStr)%spin,str(iStr)%pos,&
                                        str(jStr)%LV,str(jStr)%spin,str(jStr)%pos,&
                                          CE%eps,mapped,identical=equivalent)
      if (equivalent) then ! The structures are equivalent *and* have exactly the same concentration
         write(*,'("WARNING: *** In the input list, structures ",i4," and ",i4," appear to be the same.")') iStr, jStr
         write(*,'("WARNING:    Structure title: ",a80)') adjustl(str(iStr)%title)
         write(*,'("WARNING:    Structure title: ",a80)') adjustl(str(jStr)%title)
         write(100,'(I5,1x,I5,5x,A1,A,A,A,A1)') iStr, jStr, """", trim(adjustl(str(iStr)%title)),""" <--> """,trim(adjustl(str(jStr)%title)),""""
         !stop
      endif
   enddo
enddo
close(100)
endif ! is no Reference

! Check that each structure doesn't have two atoms on the same site
!do iStr = 1, nStr ! Loop over each structures
!  call check_for_duplicate_atoms(str(iStr),iAt,jAt,CE%eps)
!  
!  if (iAt/=0 .or. jAt/=0) then
!    if (.not. isReference) then
!      write(*,'(A,I4,A)') "ERROR: in input structure #", iStr
!    else
!      write(*,'(A,I4,A,I4,A)') "ERROR: in reference #",iRef," structure #",istr
!    endif
!    write(*,'("ERROR:   Structure has duplicate atoms on one site.")') 
!    write(*,'("ERROR:   Duplicate atoms are numbers ",i3," and ",i3)') iAt, jAt
!    write(*,'("ERROR:   Structure title: ",a80)') adjustl(str(iStr)%title)
!    stop
!  endif
!
!enddo

! Last check is to see if each structure has the correct volume/atom. If it does, and 
! it passed the checks above then all the atomic sites in the cell must be occupied.
pVol = abs(determinant(CE%pLV)) ! Volume of the parent cell
do iStr = 1, nStr ! Loop over each structure
   ! Check the volume/atom ratio of each structures with the volume of the parent cell
   ! This will not work as coded if the parent cell has more than one site/cell
!   if (.not.(equal(nPar*abs(determinant(str(iStr)%LV))/size(str(iStr)%spin,1),pVol,CE%eps))) then
!      write(*,'("There is a problem with the volume/atom ratio input structure ",i4)') iStr
!      write(*,'("Structure title: ",a80)') adjustl(str(iStr)%title)
!      stop
!   endif
enddo

! Last check. Is the number of types in the input structure consistent with the CE rank from lat.in?
! This routine won't win any elegance awards but it seems to do the trick.
allocate(uqlist(CE%CErank))
uqlist = -99  ! set it to something very unlikely
do iStr = 1, nStr ! Loop over each structure
   nAt = size(str(iStr)%pos,2)
   do iAt = 1, nAt ! Loop over each atom in the current structure
      found = .false.
      do iuq = 1, count(uqlist/=-99) ! Loop over every known (so far) atom type
         if (iuq>CE%CErank) then
            write(*,'("More atom types in input structures than in lat.in file")')
            write(*,'("Structure title: ",a80)') adjustl(str(iStr)%title)
            stop
         endif
!         write(*,'("iuq",i3,"  rank",i3,"  str",i3,"  iAt",i3," uqlist",i3,"  spin",i3)') &
!                       iuq,CE%CErank,istr,iAt,uqlist(iuq),str(iStr)%spin(iAt)
         if (uqlist(iuq)==str(iStr)%spin(iAt)) found = .true.
      enddo ! End loop over unique types found so far
      if (.not. found) uqlist(iuq) = str(iStr)%spin(iAt)
   enddo
enddo
! I don't think there are any other possibilities of bad input if all these tests pass...
ENDSUBROUTINE check_input_structures

end program convert_structures_to_enumformat
