MODULE enumeration_utilities
use num_types
use utilities
use symmetry
use vector_matrix_utilities
use numerical_utilities
use rational_mathematics
use enumeration_types
use derivative_structure_generator
use labeling_related
implicit none
private
public  map_enumStr_to_real_space, cartesian2direct, read_poscar, &
        get_HNF_of_derivative_structure,get_HNF_of_derivative_structure_old, &
        compare_two_gstructures, &
        get_gspace_representation, find_match_in_structenumout, find_equivalent_labelings

CONTAINS

  !!<summary>This subroutine maps an enumerated structure to real space.</summary>
  !!<parameter name="k" regular="true">Number of atom types (number of colors in the
  !!enumeration)</parameter>
  !!<parameter name="n" regular="tue">Number of atoms in the unit cell of given structure
  !!</parameter>
  !!<parameter name="HNF" regular="true">Hermite normal form
  !!corresponding to the superlattice </parameter>
  !!<parameter name="labeling" regular="true">List, 0..k-1, of the atomic type at
  !!each site</parameter>
  !!<parameter name="pLV" regular="true">parent lattice vectors (Cartesian coordinates)
  !!</parameter>
  !!<parameter name="pBas" regular="true">3xn array of the interior
  !!points of the parent lattice </parameter>
  !!<parameter name="eps" regular="true">Finite precision factor.</parameter>
  !!<parameter name="sLV" regular="true">Superlattice vectors (Cartesian coordinates)
  !!</parameter>
  !!<parameter name="aBas">Atomic positions (cartesian coordinates)</parameter>
  !!<parameter name="spin">Occupation variable of the positions</parameter>
  !!<parameter name="gIndx">Label ordinal for each atom (position in the labeling)</parameter>
  !!<parameter name="x">Concentration of each component</parameter>
  !!<parameter name="L" regular="true">HNF-&gt;SNF left
  !!transform. Need this to map r-&gt;G (Eq. 3 in first
  !!paper)</parameter>
  !!<parameter name="S" regular="true">Diagonal entries of the SNF</parameter>
  !!<parameter name="minkowskiReduce" regular="true">flag whether to
  !!perform a minkowski reduction or not.  (inserted by tk, need it in
  !!UNCLE). Problem is that the minkowski reduction may exchange the
  !!LVs of the supercell, which causes trouble for surface-slab
  !!reference structuresit.</parameter>
  SUBROUTINE map_enumStr_to_real_space(k, n, HNF, labeling, pLV, pBas, eps, sLV, aBas, spin, gIndx, x, L, S, minkowskiReduce)
    integer, intent(in) :: k 
    integer, intent(in) :: n 
    integer, intent(in) :: HNF(3,3) 
    character(500), intent(in) :: labeling 
    real(dp), intent(in) :: pLV(3,3), eps 
    real(dp), intent(in) :: pBas(:,:) 
    real(dp), intent(out) :: sLV(3,3) 
    real(dp), pointer :: aBas(:,:) 
    integer, pointer :: spin(:) 
    integer, pointer :: gIndx(:) 
    real(dp), pointer :: x(:) 
    integer, intent(in):: L(3,3) 
    integer, intent(in):: S(3) 
    logical, intent(in) :: minkowskiReduce 

    integer a, b, c, d, e, f ! elements of the HNF matrix
    integer digit ! one of the labelings in the labeling
    integer ic, z1, z2, z3 ! Counter over number of interior points, steps in the g table
    !("lattice coords")
    integer iAt, i, iD, nD
    real(dp) :: greal(3)   ! Floating point representation of the group element components
    !(g-vector)
    integer  :: g(3)       ! Integer version of greal

    nD = size(pBas,2)
    !allocate(gotAtomFromLabPos(n*nD)); gotAtomFromLabPos = .false.

    ! Define the non-zero elements of the HNF matrix
    a = HNF(1,1); b = HNF(2,1); c = HNF(2,2)
    d = HNF(3,1); e = HNF(3,2); f = HNF(3,3)
    ! Compute the superlattice vectors 
    !write(*,'("Parent lattice:",/,3(3f7.3,1x,/))') (pLV(:,i),i=1,3) 
    !write(*,'("Parent lattice:",/,3(3I,1x,/))') (HNF(i,:),i=1,3) 

    sLV = matmul(pLV,HNF)
    !write(*,'("Superlattice:",/,3(3f7.3,1x,/))') (sLV(:,i),i=1,3) 

    ! Find the coordinates of the basis atoms
    allocate(aBas(3,n*nD))
    allocate(spin(n*nD),gIndx(n*nD))
    gIndx=-1

    !write(*,'(3(f7.3,1x))') (sLV(i,:),i=1,3)
    ! Let's get the fattest basis (Minkowski reduction)
    if (minkowskiReduce) call minkowski_reduce_basis(sLV,sLV,eps)
    !write(*,'(3(f7.3,1x))') (sLV(i,:),i=1,3)


    ! Find each atomic position from the g-space information
    ic = 0  ! Keep track of the number of points mapped so far
    do iD = 1, nD ! Loop over the number at sites/parent cell (the d set)
       do z1 = 0, a-1 ! For the limits on the loops, see the "interior_points.pdf" write-up
          do z2 = (b*z1)/a, c+(b*z1)/a - 1
             do z3 = z1*(d-(e*b)/c)/a+(e*z2)/c, f+z1*(d-(e*b)/c)/a+(e*z2)/c - 1
                ic = ic + 1 !; if (ic > n) stop "Problem in basis atoms..."
                ! Atomic basis vector in Cartesian coordinates
                aBas(:,ic)=matmul(pLV,(/z1,z2,z3/))+pBas(:,iD)
                greal = matmul(L,(/z1,z2,z3/)) ! Map position into the group
                g = nint(greal) ! Convert the g-vector from real to integer
                if(.not. equal(greal,g,eps)) stop "map2G didn't work in map_enumStr_to_real_space"
                g = modulo(g,S) ! Bring the g-vector back into the first tile
                ! gIndx is the index in the configuration string that tells us which atom type
                !is used at this position
                gIndx(ic) = (iD-1)*S(1)*S(2)*S(3)+g(1)*S(2)*S(3)+g(2)*S(3)+g(3)+1
             enddo
          enddo
       enddo
    enddo
    if (ic /= n*nD) stop "ERROR: map_enumStr_to_real_space: Didn't find the correct # of basis atoms"


    allocate(x(k))
    x = 0.0
    if (mod(k,2)==0) then
       do iAt = 1, n*nD
          i = ichar(labeling(gIndx(iAt):gIndx(iAt)))-48
          !      gotAtomFromLabPos(gIndx(iAt)) = .true.
          digit = i-k/2 ! convert 0..k-1 label to spin variable -k/2..k/2
          x(i+1) = x(i+1) + 1  ! Keep track of the concentration of each atom type
          if (digit<0) then
             spin(iAt) = digit
          else
             spin(iAt) = digit+1 ! skip 0 as a spin if k is even
          endif
       enddo
    else
       do iAt = 1, n*nD
          i = ichar(labeling(gIndx(iAt):gIndx(iAt)))-48
          spin(iAt) = i-k/2
          x(i+1) = x(i+1) + 1 ! Keep track of the concentration of each atom type
       enddo
    endif
    x = x/real(n*nD,dp)

    !if (.not. all(gotAtomFromLabPos .eqv. .true.)) stop "Labeling to atom conversion failed
    !in map_enumStr_to_real_space"
    !deallocate(gotAtomFromLabPos)
  ENDSUBROUTINE map_enumStr_to_real_space

  !!<summary>This routine takes three lattice vectors and a list of
  !!atomic basis vector in Cartesian coordinates and converts them to
  !!direct ("lattice") coordinates.</summary>
  !!<parameter name="sLV" regular="true">Superlattice vectors
  !!(Cartesian coordinates)</parameter>
  !!<parameter name="aBas" regular="true">Atomic positions (cartesian
  !!coordinates first, then direct)</parameter>
  !!<parameter name="eps" regular="true">Finite precision
  !!tolerance.</parameter>
  subroutine cartesian2direct(sLV,aBas, eps) 
    real(dp), intent(in)    :: sLV(3,3)  
    real(dp), intent(inout) :: aBas(:,:) 
    real(dp), intent(in)    :: eps

    real(dp) :: sLVinv(3,3)
    integer iAt, nAt

    nAt = size(aBas,2)
    call matrix_inverse(sLV,sLVinv)
    ! Convert aBas to DIRECT COORDINATES
    do iAt=1,nAt
       aBas(:,iAt) = matmul(sLVinv,aBas(:,iAt)) ! Put positions into
       ! "direct" coordinates This keeps the atomic coordinates inside
       ! the first unit cell--- not necessary but aesthetically
       ! pleasing.
       do while(any(aBas(:,iAt) >= 1.0_dp - eps) .or. any(aBas(:,iAt) < 0.0_dp - eps)) 
          aBas(:,iAt) = merge(aBas(:,iAt), aBas(:,iAt) - 1.0_dp, aBas(:,iAt) <  1.0_dp - eps) 
          aBas(:,iAt) = merge(aBas(:,iAt), aBas(:,iAt) + 1.0_dp, aBas(:,iAt) >= 0.0_dp - eps)
       enddo
    end do
    ! End of conversion to DIRECT COORDINATES
  end subroutine cartesian2direct

  !!<summary>Reads data in from a poscar file.</summary>
  !!<parameter name="fname" regular="true">Name of file to be read.</parameter>
  !!<parameter name="title" regular="true"></parameter>
  !!<parameter name="sLV" regular="true">Superlattice vector of the
  !!read-in structure</parameter>
  !!<parameter name="aBas">The atomic basis vectors of the atoms (in
  !!Cartesian, not direct coords)</parameter>
  !!<parameter name="aTyp">the atomic label (ie, type) for each atom
  !!in the atomic basis.</parameter>
  SUBROUTINE read_poscar(fname,title,sLV,aBas,aTyp)
    character(80), intent(in) :: fname
    character(80), intent(out):: title
    real(dp), intent(out) :: sLV(3,3) 
    real(dp), pointer :: aBas(:,:) ! (out) 
    integer, pointer :: aTyp(:) ! (out) 
    
    integer iAt, iLV, status, ic, nTyp, nAt, idx
    real(dp) scale, sLVinv(3,3)
    character(80) coords, dummy, d2 ! Direct or Cartesian; input_string
    logical err
    integer, allocatable :: aList(:)
    
    open(16,file=fname,iostat=status)
    if(status/=0)then; write(*,'("The file specified for read_poscar in enumeration_utilities does not exist")')
       write(*,'("Filename: ",a80)') fname; stop; endif
    open(17,file="readcheck_poscar_check.out")
    
    read(16,'(a80)') title
    write(17,'(a80)') title
    read(16,*) scale
    write(17,'("Scale factor: ",f8.5,/)') scale
    read(16,*) (sLV(:,iLV),iLV=1,3)
    write(17,'("Vectors (as columns)",/,3(3(f7.3,1x),/))') (sLV(iLV,:),iLV=1,3)
    sLV = sLV*scale
    write(17,'("Vectors (scaled)",/,3(3(f7.3,1x),/))') (sLV(iLV,:),iLV=1,3)
    call matrix_inverse(sLV,sLVinv,err)
    if (err) stop "Co-planar vectors in call to read_poscar"
    read(16,'(a80)') dummy ! list of # of each atom type
       
    ! Need to parse this line into a bunch of numbers and store them in an array
    dummy=trim(adjustl(dummy)) ! Get rid of leading and trailing blanks
    d2 = dummy
    ! Now find out how many number there are in the dummy line
    ic = 0
    do
       ic = ic + 1
       read(d2,*) status
       d2=d2(index(d2," "):80)
       d2 = adjustl(d2)
       if(len_trim(d2)==0) exit
       if(ic>25) stop "too many atoms in POSCAR: read_poscar in enumeration_utilities" 
    enddo
    nTyp = ic
    allocate(aList(nTyp))
    write(17,'(" Number of each kind of atom:")')
    do ic = 1, nTyp
       read(dummy,*) aList(ic)
       dummy = dummy(index(dummy," "):80)
       write(17,'(i3)',advance="no") aList(ic)
    enddo
    write(17,*)
    nAt = sum(aList)
    allocate(aBas(3,nAt),aTyp(nAt))
    write(17,'("Number of atoms in structure: ",i3)') nAt
       
    read(16,*) coords
    call ucase(coords)
    write(17,'("Coordinates (direct or Cartesian): ",a1)') coords
    ! Read in all the atoms in the POSCAR. Convert to Cartesian coordinates, if necessary
    do iAt = 1, nAt
       read(16,*) aBas(:,iAt)
       write(17,'("Atom #: ",i3," position: ",3(f7.3,1x))') iAt,aBas(:,iAt) 
       if(coords(1:1) .EQ. "D") then ! convert to Cartesian
          aBas(:,iAt) = matmul(sLV,aBas(:,iAt))
          write(17,'("Atom #: ",i3," position (Cart): ",3(f7.3,1x))') iAt,aBas(:,iAt)
       endif
    enddo
    close(16)
    ! Set the type for each atom
    idx = 1
    do iAt = 1, nAt
       aTyp(iAt) = idx-1
       write(17,'("#:",i3,"   idx:",i3)') iAt, idx
       if(iAt == sum(aList(1:idx))) idx = idx + 1
    enddo
    close(17)
  END SUBROUTINE read_poscar

  !!<summary>This subroutine takes a structure defined in real space
  !!(atomic basis vectors in Cartesian coordinates) and checks whether
  !!or not it is a derivative structure of a parent lattice. If it is,
  !!the parent lattice, d-set, and HNF/SNF are returned.</summary>
  !!<parameter name="title" regular="true"></parameter>
  !!<parameter name="sLV" regular="true">Input lattice
  !!vectors</parameter>
  !!<parameter name="aBas" regular="true">Input atomic coordinates
  !!(only changed if a shift is needed*) </parameter>
  !!<parameter name="aTyp" regular="true">Type of each
  !!atom</parameter>
  !!<parameter name="pLV" regular="true">Lattice vectors of the parent
  !!lattice</parameter>
  !!<parameter name="dset">d-set of the parent lattice</parameter>
  !!<parameter name="HNF">List of (rotationally-) equivalent
  !!HNFs</parameter>
  !!<parameter name="SNF" regular="true"></parameter>
  !!<parameter name="L" regular="true">L is the left transform for
  !!SNF</parameter>
  !!<parameter name="eps" regular="true">Finite precision
  !!tolerance.</parameter>
  !!<parameter name="fixedOcc">List of the atoms with fixed
  !!occupation. This is important for the symmetry and to determine
  !!whether a cell is primitive</parameter>
  SUBROUTINE get_HNF_of_derivative_structure(title,sLV,aBas,aTyp,pLV,dset,HNF,SNF,L,eps,fixedOcc)
    character(80), intent(in) :: title
    real(dp), intent(in)      :: sLV(3,3)       
    real(dp), intent(inout)   :: aBas(:,:)      
    integer, intent(in)       :: aTyp(:)        
    real(dp), intent(in)      :: pLV(:,:)       
    real(dp), pointer         :: dset(:,:)      ! (intent(in)) 
    real(dp), pointer         :: dsetStruc(:,:) ! d-set of structure
    ! (is determined from the structure, should in the end be
    ! equivalent to dset)
    integer                   :: nD
    real(dp), intent(in)      :: eps
    integer, intent(in), target, optional :: fixedOcc(:)  
    integer, pointer                      :: fixedOcc_(:) ! internal variable for fixedOcc
    integer, pointer                      :: removedAtoms(:)  ! indices
    ! of those atoms that were removed by make_primitive
    integer, pointer :: HNF(:,:,:) ! (out) 
    integer, intent(out) :: SNF(3,3), L(3,3) 

    integer nAt, iAt, nOp, iOp, row(6), iuq, nuq, i, iE, iP
    integer, pointer :: aTypTemp(:)
    real(dp), pointer :: aBasTemp(:,:), sgrot(:,:,:), sgshift(:,:)
    real(dp), dimension(3,3) :: pLVinv, pLVtemp, parLattTest, sLVinv
    integer,dimension(3,3) :: S, R, tempH, newH
    logical err, unique, mapped
    integer, allocatable :: trow(:,:), vs(:), idx(:)
    real(dp) :: testsite(3), diff(3)

    nAt = size(aBas,2)
    if(nAt/=size(aTyp)) stop "Input to get_HNF_of_derivative_structure is inconsistent"

    nD = size(dset,2)
    
    if (present(fixedOcc)) then; fixedOcc_ => fixedOcc
    else; allocate(fixedOcc_(0))
    endif

    ! open(17,file="debug_get_HNF.out")
    ! write(17,'(3/,"<<< Finding the g-space representation >>>",2/)')

    ! do iAt = 1, nAt
    !    write(17,'("Atom #: ",i2," position: ",3(f7.3,1x))') iAt, aBas(:,iAt)
    ! enddo

    allocate(aTypTemp(nAt),aBasTemp(3,nAt))
    aTypTemp = aTyp; aBasTemp = aBas
    pLVtemp = sLV 
    call make_primitive(pLVtemp,aTypTemp,aBasTemp,.false.,eps,removed_=removedAtoms)
    if(nAt/=size(aTypTemp)) then; ! the structure was reduced by make_primitive
       ! write(17,'(/,"ERROR: The input structure wasn''t primitive")')
       ! write(17,'("atom type: ",20(i2,1x))') aTypTemp(:)
       ! write(17,'("number of atoms: ",i2,5x," size of aTyp:",i2)') nAt, size(aTypTemp)
       !stop "ERROR: input structure for get_HNF_of_derivative_structure was not primitive";endif
       write(*,'(/,A)') "WARNING: input structure for get_HNF_of_derivative_structure was not primitive"; 
       allocate(hnf(1,1,1)); HNF = 0
       return
    endif
    
    ! Now make every atom the same and apply make_primitive to find the parent cell
    ! (but don't change the labeling of fixed occupation sites)
    aTypTemp = 1
    aTypTemp(fixedOcc_) = aTyp(fixedOcc_)
    pLVtemp = sLV
    call make_primitive(pLVtemp,aTypTemp,aBasTemp,.false.,eps)
    ! write(17,'("Parent lattice of input superlattice (columns): ",/,3(3(f8.4,1x),/))') (pLVtemp(iAt,:),iAt=1,3)

    call matrix_inverse(pLVtemp,parLattTest,err)
    if(err) stop "Problem inverting parent lattice basis"
    if (nD<size(aTypTemp))then;
       write(*,'(/,A)') "WARNING: maybe your epsilon is too small"
       allocate(hnf(1,1,1)); HNF = 0
       return
    else if (nD>size(aTypTemp))then;
    ! Note: we should try to find a better solution here. The problem
    ! here is that we can define a lattice with >1 dvector that can be
    ! effectively reduced to a different lattice with only 1
    ! dvector. If this happens, all structures end up in this "if"
    ! statement with aTypTemp<nD.
       write(*,'(/,A)') "WARNING: parent lattice does not seem to be primitive"
       allocate(hnf(1,1,1)); HNF = 0
       return
    endif
    ! If we pass this test, the two cell bases are equivalent even if
    ! not equal. From this point on, use the one from the
    ! struct_enum.out file. (that one is pLV, not pLVtemp)


    allocate(dsetStruc(3,nD)); dsetStruc = aBasTemp
    call matrix_inverse(pLV,pLVinv,err) 
    if(err) stop "matrix inverse failed in get_HNF_of_derivative_structure"
    ! write(17,'("d-set: ",/,200(3(f8.4,1x),/))') (dsetStruc(:,iAt),iAt=1,size(dsetStruc,2))
    ! write(17,'("Size of supercell: ",i3)') abs(nint(determinant(sLV)/determinant(pLV)))

    ! I think that both of the d-sets (one from struct_enum.out and
    ! one from POSCAR) should be brought the unit cell (if they aren't
    ! already) and it should be the *same* unit cell

    ! write(17,'("After bringing into common unit cell:")')
    do iAt = 1, nD
       call bring_into_cell(dsetStruc(:,iAt),pLVinv,pLV,eps) 
       call bring_into_cell(dset     (:,iAt),pLVinv,pLV,eps) 
       ! write(17,'("  dset struc: ",3(F8.4,1x))') dsetStruc(:,iAt)
       ! write(17,'("  dset      : ",3(F8.4,1x))') dset(:,iAt)   
    end do; 
    ! write(17,*)

    ! Check that the site basis for both the enum file and POSCAR are
    !  consistent.  They might not have the same origin so try all
    !  distinct origin shifts and see if there is any shift that maps
    !  all of the sites in one case to all of the sites in the other
    ! a) always start with 0 shift first
    !    (note: Gus's original version started with part b below. For
    !    surfaces it produced incorrect shifts, though. This is
    !    probably due to another fixed occupation issue, but tk
    !    decided not to digg to deep into it for the time being. After
    !    all, the input structures that are compared have already been
    !    checked for occupying the lattice, so a shift should not be
    !    necessary in UNCLE.)
    diff = 0._dp
    do iE = 1, nD
       mapped=.false.
       do iP = 1, nD
          testsite = dset(:,iP)-diff
          call bring_into_cell(testsite,pLVinv,pLV,eps) ! Make sure we stay in the same cell
          if (equal(testsite,dset(:,iE),eps)) then ! this site maps to one in other lattice
             mapped=.true.;
             exit
          endif
       enddo
       if(.not. mapped) exit
    enddo
    ! b) if the 0 shift was not successfull, try everything:
    if (.not. mapped) then
       do iAt = 1, nD
          !write(*,'("iAt:",i3,"  EnumBas ",3(F8.4,1x))') iAt,EnumBas(:,iAt)
          !write(*,'("iAt:",i3,"  dset ",3(F8.4,1x))') iAt,dset(:,iAt)
          diff = dsetStruc(:,iAt)-dset(:,1)
          !write(*,'("iAt:",i3,"  diff ",3(F8.4,1x))') iAt,dset(:,iAt)-EnumBas(:,1)
          do iE = 1, nD
             mapped=.false.
             do iP = 1, nD
                testsite = dset(:,iP)-diff
                !write(*,'("iP:",i3,"  testsite",3(F8.4,1x))') iP,testsite
                call bring_into_cell(testsite,pLVinv,pLV,eps) ! Make
                ! sure we stay in the same cell write(*,'("iP:",i3,"
                ! testsite",3(F8.4,1x))') iP,testsite
                if (equal(testsite,dset(:,iE),eps)) then ! this site
                   ! maps to one in other lattice
                   mapped=.true.;
                   exit
                endif
             enddo
             if(.not. mapped) exit
          enddo
          ! If we get to here and mapped is true, then the iE loop
          ! concluded without mapped ever coming up false in the iP
          ! loop. If that is the case, then the sites in EnumBas all had
          ! coincident sites in aBas. That means that the parent lattices
          ! (POSCAR and enum) are equivalent even if they don't have the
          ! same origin.
          if (mapped) exit
       enddo
    endif

    if (.not. mapped) stop "The lattice sites of the two parent lattices are not coincident"
    if (.not. equal(diff,0._dp,eps)) then
       ! write(17,'("The atomic sites of the parent lattice of the input structure have been shifted by")')
       ! write(17,'(3(f8.4,1x),/)') diff
       ! write(*,'("DEBUG: The atomic sites of the parent lattice of the input structure have been shifted by")')
       ! write(*,'("DEBUG: ",3(f8.4,1x),/)') diff

       !Gus>   dset = dset - spread(diff,2,size(aBas(:,2)))
       dset = dset - spread(diff,2,size(aBas,2))

       do iAt = 1, nD
          call bring_into_cell(dset(:,iAt),pLVinv,pLV,eps)
          ! write(17,'("new pBas: ",3(f8.4,1x))') dset(:,iAt)
       enddo
       ! write(17,'(/)')
       call matrix_inverse(sLV,sLVinv,err)
       if(err) stop "Superlattice vectors are co-planar in routine ""get_HNF_of_derivative_superstructure"""
       do iAt = 1, size(aBas,2)
          aBas(:,iAt) = aBas(:,iAt) - diff
          call bring_into_cell(aBas(:,iAt),sLVinv,sLV,eps)
          ! write(17,'("new sLV bas: ",i2,3x,3(f8.4,1x))') iAt, aBas(:,iAt)
       enddo
       ! write(17,'(//)')
    endif
    
    if(.not. equal(determinant(sLV)/determinant(pLV),nint(determinant(sLV)/determinant(pLV)),eps)) stop &
         "The superlattice in not an integer multiple of the primitive"
    if(.not. equal(matmul(pLVinv,sLV),nint(matmul(pLVinv,sLV)),eps)) stop &
         "ERROR: HNF was non-integer"
    S = nint(matmul(pLVinv,sLV))

    call get_spaceGroup(pLV,aTypTemp,dset,sgrot,sgshift,.false.,eps)
    ! I think it's OK here not to pass in labels for the types
    ! (aTypTemp=1). If we pick up extra symmetries from d-set
    ! permutations, these will not be effective if the allowed labels
    ! are different. In other words, it can't hurt anything.

    ! write(17,'("After get_spaceGroup:")')
    do iAt = 1, nD
       call bring_into_cell(dsetStruc(:,iAt),pLVinv,pLV,eps) 
       call bring_into_cell(dset     (:,iAt),pLVinv,pLV,eps) 
       ! write(17,'("  dset struc: ",3(F8.4,1x))') dsetStruc(:,iAt)
       ! write(17,'("  dset      : ",3(F8.4,1x))') dset(:,iAt)   
    end do;

    nOp = size(sgrot,3)
    allocate(trow(nOp,6),vs(nOp),idx(nOp))
    ! open(16,file="debug_get_spacegroup.out")
    ! do iOp = 1, nOp
    !    do i = 1, 3
    !       write(16,'("Op#: ",i2,3x,3(f8.4,1x))') iOp,sgrot(i,:,iOp)
    !    end do
    !    write(16,'("Shift#: ",3(f8.4,1x),/)') sgshift(:,iOp)
    ! end do

    ! Find the unrotated form of this structure's HNF. We'll need this
    ! later to get the SNF's left transform matrix, L, to map atom
    ! postions into the group. The find_gspace_representation routine
    ! will assume that the *first* HNF in the list in the unrotated
    ! one.
    call HermiteNormalForm(S,newH,R)
    trow(1,:) = (/newH(1,1),newH(2,1),newH(2,2),newH(3,1),newH(3,2),newH(3,3)/)
    
    idx = (/(i,i=1,nOp)/)
    nuq = 1
    do iOp = 1, nOp
       ! Under the current rotation, what HNF does the current HNF turn into?
       tempH = nint(matmul(pLVinv,matmul(sgrot(:,:,iOp),matmul(pLV,S))))
       call HermiteNormalForm(tempH,newH,R)
       unique = .true.
       row = (/newH(1,1),newH(2,1),newH(2,2),newH(3,1),newH(3,2),newH(3,3)/)
       ! write(17,'("HNF (rot #):",i3,3x,6(i2,1x))') iOp,row
       do iuq = 1, nuq
          if (all(row == trow(iuq,:))) then
             unique = .false.
             exit
          endif
       enddo
       if (unique) then
          nuq = nuq + 1
          trow(nuq,:) = row
       endif
    enddo
    ! write(17,'("Unique HNF entries after rotations:")')
    ! do i = 1,nuq
    !    write(17,'("row #:",i3,3x,6(i2,1x))') i,trow(i,:)
    ! enddo
    allocate(HNF(3,3,nuq))
    HNF = 0
    do iuq = 1,nuq
       HNF(1,1,iuq) = trow(iuq,1); HNF(2,1,iuq) = trow(iuq,2); HNF(2,2,iuq) = trow(iuq,3)
       HNF(3,1,iuq) = trow(iuq,4); HNF(3,2,iuq) = trow(iuq,5); HNF(3,3,iuq) = trow(iuq,6)
       ! write(17,'("HNF #: ",i3,/,3(3(i2,1x),/))') iuq,(HNF(iAt,:,iuq),iAt=1,3)
    enddo
    
    call SmithNormalForm(S,L,SNF,R)
    ! write(17,'("Integer transform of superlattice: ",/,3(3(i2,1x),/))') (S(iAt,:),iAt=1,3)
    ! write(17,'("SNF: ",/,3(3(i2,1x),/))') (SNF(iAt,:),iAt=1,3)
    ! close(17)
    
  END SUBROUTINE get_HNF_of_derivative_structure

  !!<summary>This subroutine takes a structure defined in real space
  !!(atomic basis vectors in Cartesian coordinates) and checks whether
  !!or not it is a derivative structure of a parent lattice. If it is,
  !!the parent lattice, d-set, and HNF/SNF are returned.
  !!<comments>THIS IS GUS'S ORIGINAL. tk HAS MODIFIED IT IN ORDER TO
  !!GET THE STRUCTURE COMPARISON THIS ROUTINE IS ONLY NEEDED IN THE
  !!DRIVER ROUTINE FIND_STRUCTURE_IN_LIST.X, and this modified version
  !!is in GET_HNF_OF_DERIATIVE_STRUCTURE above. PLEASE: Don't change
  !!anything here that could as well be of importance to
  !!get_hnf_of_derivative_structure above.</comments></summary>
  !!<parameter name="sfname" regular="true">Name of file to be
  !!searched for target structure.</parameter>
  !!<parameter name="sLV" regular="true">Input lattice
  !!vectors.</parameter>
  !!<parameter name="aBas" regular="true">Input atomic coordinates
  !!(only changed if a shift is needed*).</parameter>
  !!<parameter name="aTyp" regular="true">Type of each atom.</parameter>
  !!<parameter name="pLV" regular="true">Lattice vectors of the parent
  !!lattice.</parameter>
  !!<parameter name="dset">d-set of the parent lattice</parameter>
  !!<parameter name="HNF">List of (rotationally-) equivalent
  !!HNFs</parameter>
  !!<parameter name="SNF" regular="true"></parameter>
  !!<parameter name="L" regular="true">L is the left transform for
  !!SNF.</parameter>
  !!<parameter name="eps" regular="true">Finite precision tolerance.</parameter>
  SUBROUTINE get_HNF_of_derivative_structure_old(sfname,sLV,aBas,aTyp,pLV,dset,HNF,SNF,L,eps)
    character(80), intent(in) :: sfname 
    real(dp), intent(in) :: sLV(3,3)    
    real(dp), intent(inout):: aBas(:,:) 
    integer, intent(in) :: aTyp(:) 
    real(dp), intent(out):: pLV(:,:) 
    real(dp), pointer :: dset(:,:) ! (out)  
    real(dp), intent(in):: eps
    integer, pointer :: HNF(:,:,:) ! (out) 
    integer, intent(out) :: SNF(3,3), L(3,3) 

    integer nAt, iAt, nOp, iOp, row(6), iuq, nuq, i, nEnumBas, iE, iP, ioerr
    integer, pointer :: aTypTemp(:)
    real(dp), pointer :: aBasTemp(:,:), sgrot(:,:,:), sgshift(:,:), EnumBas(:,:)
    real(dp), dimension(3,3) :: pLVinv, pLVtemp, parLattTest, sLVinv
    integer,dimension(3,3) :: S, R, tempH, newH
    logical err, unique, mapped
    integer, allocatable :: trow(:,:), vs(:), idx(:)
    real(dp) :: testsite(3), diff(3)

    nAt = size(aBas,2)
    if(nAt/=size(aTyp)) stop "Input to get_HNF_of_derivative_structure is inconsistent"
    
    open(17,file="debug_get_HNF.out")
    write(17,'(3/,"<<< Finding the g-space representation >>>",2/)')
    
    do iAt = 1, nAt
       write(17,'("Atom #: ",i2," position: ",3(f7.3,1x))') iAt, aBas(:,iAt)
    enddo
    
    ! Need to read in the parent lattice from struct_enum.out. The
    ! parent lattice of struct_enum.out and the parent lattice of the
    ! structure to be checked must, of course, be equivalent. But for
    ! checking the structures via the HNF, the same parent lattice
    ! representation must be used for both. (That is, the bases for
    ! the parent lattice must not only be equivalent, but identical.)
    ! This is because, if A_1!=A_2 (the two bases for the parent
    ! lattice), then A_1*H!=A_2*H. We need to use the same parent
    ! basis when extracting the HNF of a superlattice.
    open(18,file=sfname,status="old",iostat=ioerr)
    read(18,*); read(18,*) ! Skip the first two lines
    do i = 1, 3; read(18,*) pLV(:,i); enddo
    write(17,'("parent lattice vectors (columns) from ",a80)') adjustl(sfname)
    write(17,'(3(3f8.4,1x,/))') transpose(pLV)

    ! Need to make sure that the interior points of the parent
    ! lattice are also equivalent so read them in as well.
    read(18,*) nEnumBas
    allocate(EnumBas(3,nEnumBas))
    write(17,'("Interior points of the parent (multi)lattice: ")')
    do iAt = 1, nEnumBas
       read(18,*) EnumBas(:,iAt)
       write(17,'("pBas: ",3(f8.4,1x))') EnumBas(:,iAt)
    enddo
       
    allocate(aTypTemp(nAt),aBasTemp(3,nAt))
    aTypTemp = aTyp; aBasTemp = aBas
    pLVtemp = sLV 
    call make_primitive(pLVtemp,aTypTemp,aBasTemp,.false.,eps)
    if(nAt/=size(aTypTemp)) then; ! the structure was reduced by make_primitive
       write(17,'(/,"ERROR: The input structure wasn''t primitive")')
       write(17,'("atom type: ",20(i2,1x))') aTypTemp(:)
       write(17,'("number of atoms: ",i2,5x," size of aTyp:",i2)') nAt, size(aTypTemp)
       stop "ERROR: input structure for get_HNF_of_derivative_structure was not primitive"
    endif
    aTypTemp = 1; pLVtemp = sLV
    ! Now make every atom the same and apply make_primitive to find the parent cell
    call make_primitive(pLVtemp,aTypTemp,aBasTemp,.false.,eps)
    !write(17,'("Parent lattice of input superlattice (columns): ",/,3(3(f8.4,1x),/))') (pLVtemp(iAt,:),iAt=1,3)
    call matrix_inverse(pLVtemp,parLattTest,err)
    if(err) stop "Problem inverting parent lattice basis"
    if (nEnumBas/=size(aTypTemp))then;print*,"Number of parent lattice sites in "//trim(adjustl(sfname))
       print *,"isn't the same as the structure being checked";
       print *,"size(aTypTemp)",size(aTypTemp),"nEnumBas",nEnumBas
       print *,"Maybe your input cell did not reduce to the same-sized primitive as in the struct_enum.out"
       print *,"routine: get_HNF_of_derivative_structure_old, file: enumeration_utilities.f90";stop
    endif
    ! If we pass this test, the two cell bases are equivalent even
    ! if not equal. From this point on, use the one from the
    ! struct_enum.out file. (that one is pLV, not pLVtemp)

    allocate(dset(3,size(aTypTemp))); dset = aBasTemp
    call matrix_inverse(pLV,pLVinv,err) 
    if(err) stop "matrix inverse failed in get_HNF_of_derivative_structure"
    write(17,'("d-set: ",/,200(3(f8.4,1x),/))') (dset(:,iAt),iAt=1,size(dset,2))
    write(17,'("Size of supercell: ",i3)') abs(nint(determinant(sLV)/determinant(pLV)))
       
    ! I think that both of the d-sets (one from struct_enum.out and
    ! one from POSCAR) should be brought the unit cell (if they
    ! aren't already) and it should be the *same* unit cell
    write(17,'("After bringing into common unit cell:")')
    do iAt = 1, nEnumBas
       call bring_into_cell(dset(:,iAt),pLVinv,pLV,eps) 
       call bring_into_cell(EnumBas(:,iAt),pLVinv,pLV,eps) 
       write(17,'("  EnumBas: ",3(F8.4,1x))') EnumBas(:,iAt)
       write(17,'("     dset: ",3(F8.4,1x))') dset(:,iAt)   
    end do; write(17,*)
       
    !  Check that the site basis for both the enum file and POSCAR
    !  are consistent.  They might not have the same origin so try
    !  all distinct origin shifts and see if there is any shift
    !  that maps all of the sites in one case to all of the sites
    !  in the other
    do iAt = 1, nEnumBas
       !write(*,'("iAt:",i3,"  EnumBas ",3(F8.4,1x))') iAt,EnumBas(:,iAt)
       !write(*,'("iAt:",i3,"  dset ",3(F8.4,1x))') iAt,dset(:,iAt)
       diff = dset(:,iAt)-EnumBas(:,1)
       !write(*,'("iAt:",i3,"  diff ",3(F8.4,1x))') iAt,dset(:,iAt)-EnumBas(:,1)
       do iE = 1, nEnumBas
          mapped=.false.
          do iP = 1, nEnumBas
             testsite = dset(:,iP)-diff
             !write(*,'("iP:",i3,"  testsite",3(F8.4,1x))') iP,testsite
             call bring_into_cell(testsite,pLVinv,pLV,eps) ! Make sure we
             ! stay in the same cell        
             !write(*,'("iP:",i3,"  testsite",3(F8.4,1x))') iP,testsite
             if (equal(testsite,EnumBas(:,iE),eps)) then ! this site maps
                ! to one in other lattice
                mapped=.true.;
                exit
             endif
          enddo
          if(.not. mapped) exit
       enddo
       ! If we get to here and mapped is true, then the iE loop
       ! concluded without mapped ever coming up false in the iP
       ! loop. If that is the case, then the sites in EnumBas all
       ! had coincident sites in aBas. That means that the parent
       ! lattices (POSCAR and enum) are equivalent even if they
       ! don't have the same origin.
       if (mapped) exit
    enddo
    if (.not. mapped) stop "The lattice sites of the two parent lattices are not coincident"
    if (.not. equal(diff,0._dp,eps)) then
       write(17,'("The atomic sites of the parent lattice of the input structure have been shifted by")')
       write(17,'(3(f8.4,1x),/)') diff
       dset = dset - spread(diff,2,size(aBas(:,2)))
       do iAt = 1, nEnumBas
          call bring_into_cell(dset(:,iAt),pLVinv,pLV,eps)
          write(17,'("new pBas: ",3(f8.4,1x))') dset(:,iAt)
       enddo
       write(17,'(/)')
       call matrix_inverse(sLV,sLVinv,err)
       if(err) stop "Superlattice vectors are co-planar in routine ""get_HNF_of_derivative_superstructure"""
       do iAt = 1, size(aBas,2)
          aBas(:,iAt) = aBas(:,iAt) - diff
          call bring_into_cell(aBas(:,iAt),sLVinv,sLV,eps)
          write(17,'("new sLV bas: ",i2,3x,3(f8.4,1x))') iAt, aBas(:,iAt)
       enddo
       write(17,'(//)')
    endif
       
    if(.not. equal(determinant(sLV)/determinant(pLV),nint(determinant(sLV)/determinant(pLV)),eps)) stop &
         "The superlattice in not an integer multiple of the primitive"
    if(.not. equal(matmul(pLVinv,sLV),nint(matmul(pLVinv,sLV)),eps)) stop &
         "ERROR: HNF was non-integer"
    S = nint(matmul(pLVinv,sLV))
    
    call get_spaceGroup(pLV,aTypTemp,dset,sgrot,sgshift,.false.,eps)
    ! I think it's OK here not to pass in labels for the types
    ! (aTypTemp=1). If we pick up extra symmetries from d-set
    ! permutations, these will not be effective if the allowed
    ! labels are different. In other words, it can't hurt anything.
    
    write(17,'("After get_spaceGroup:")')
    do iAt = 1, nEnumBas
       call bring_into_cell(dset(:,iAt),pLVinv,pLV,eps) 
       call bring_into_cell(EnumBas(:,iAt),pLVinv,pLV,eps) 
       write(17,'("  EnumBas: ",3(F8.4,1x))') EnumBas(:,iAt)
       write(17,'("     dset: ",3(F8.4,1x))') dset(:,iAt)   
    end do;
    
    nOp = size(sgrot,3)
    allocate(trow(nOp,6),vs(nOp),idx(nOp))
    open(16,file="debug_get_spacegroup.out")
    do iOp = 1, nOp
       do i = 1, 3
          write(16,'("Op#: ",i2,3x,3(f8.4,1x))') iOp,sgrot(i,:,iOp)
       end do
       write(16,'("Shift#: ",3(f8.4,1x),/)') sgshift(:,iOp)
    end do

    ! Find the unrotated form of this structure's HNF. We'll need
    ! this later to get the SNF's left transform matrix, L, to map
    ! atom postions into the group. The find_gspace_representation
    ! routine will assume that the *first* HNF in the list in the
    ! unrotated one.
    call HermiteNormalForm(S,newH,R)
    trow(1,:) = (/newH(1,1),newH(2,1),newH(2,2),newH(3,1),newH(3,2),newH(3,3)/)
    
    idx = (/(i,i=1,nOp)/)
    nuq = 1
    do iOp = 1, nOp
       ! Under the current rotation, what HNF does the current HNF turn into?
       tempH = nint(matmul(pLVinv,matmul(sgrot(:,:,iOp),matmul(pLV,S))))
       call HermiteNormalForm(tempH,newH,R)
       unique = .true.
       row = (/newH(1,1),newH(2,1),newH(2,2),newH(3,1),newH(3,2),newH(3,3)/)
       write(17,'("HNF (rot #):",i3,3x,6(i2,1x))') iOp,row
       do iuq = 1, nuq
          if (all(row == trow(iuq,:))) then
             unique = .false.
             exit
          endif
       enddo
       if (unique) then
          nuq = nuq + 1
          trow(nuq,:) = row
       endif
    enddo
    write(17,'("Unique HNF entries after rotations:")')
    do i = 1,nuq
       write(17,'("row #:",i3,3x,6(i2,1x))') i,trow(i,:)
    enddo
    allocate(HNF(3,3,nuq))
    HNF = 0
    do iuq = 1,nuq
       HNF(1,1,iuq) = trow(iuq,1); HNF(2,1,iuq) = trow(iuq,2); HNF(2,2,iuq) = trow(iuq,3)
       HNF(3,1,iuq) = trow(iuq,4); HNF(3,2,iuq) = trow(iuq,5); HNF(3,3,iuq) = trow(iuq,6)
       write(17,'("HNF #: ",i3,/,3(3(i2,1x),/))') iuq,(HNF(iAt,:,iuq),iAt=1,3)
    enddo
    
    call SmithNormalForm(S,L,SNF,R)
    write(17,'("Integer transform of superlattice: ",/,3(3(i2,1x),/))') (S(iAt,:),iAt=1,3)
    write(17,'("SNF: ",/,3(3(i2,1x),/))') (SNF(iAt,:),iAt=1,3)
    close(17)
  END SUBROUTINE get_HNF_of_derivative_structure_old
  
  !!<summary>Maps a list of real space atomic basis vectors and their
  !!labels into g-space and extracts the labeling</summary>
  !!<parameter name="L" regular="true"></parameter>
  !!<parameter name="A" regular="true"></parameter>
  !!<parameter name="aBas"></parameter>
  !!<parameter name="dset" regular="true"></parameter>
  !!<parameter name="aTyp" regular="true"></parameter>
  !!<parameter name="SNF" regular="true"></parameter>
  !!<parameter name="eps" regular="true">Finite precision tolerance.</parameter>
  !!<parameter name="labeling"></parameter>
  SUBROUTINE find_labeling_from_atom_basis(L,A,aBas,dset,aTyp,SNF,eps,labeling)
    integer, dimension(3,3), intent(in) :: L
    real(dp), dimension(3,3),intent(in) :: A
    real(dp), dimension(:,:), pointer   :: aBas ! (in)
    real(dp), dimension(:,:),intent(in) :: dset
    integer, dimension(:),   intent(in) :: aTyp
    integer, dimension(3,3), intent(in) :: SNF
    real(dp),                intent(in) :: eps 
    integer, pointer                    :: labeling(:) ! out

    integer, dimension(3,size(aBas,2)) :: g
    integer, pointer :: p(:,:), p_with_d(:,:)
    integer :: diag(3), perm(size(aBas,2)), dmember(size(aBas,2))
    integer :: iAt, nAt, iD, nD, pLindx, n
    real(dp) :: Ainv(3,3), T(3,3), gtemp(3)
    logical err, mapped

    !open(18,file="debug_find_labeling.out")
    !write(18,'(3/,"Takes a list of vectors, maps them into the group, determines the labeling",2/)')

    call matrix_inverse(A,Ainv,err)
    if(err) stop "Coplanar vectors in find_labeling_from_atom_basis"
    
    diag(1) = SNF(1,1); diag(2) = SNF(2,2); diag(3) = SNF(3,3) 
    nAt = size(aBas,2)
    nD = size(dset,2) ! Number of d-vectors in the dset
    
    !write(18,'("Diagonal entries of HNF: ",3(i3,1x))') diag
    !write(18,'("Left transform:",/,3(3i3,1x,/),/)') transpose(L)
    !write(18,'("parent lattice vectors (columns):",/,3(3(f7.3,1x),/),/)') transpose(A)
    !write(18,'("parLatt inverse (columns):",/,3(3(f7.3,1x),/),/)') transpose(Ainv)
    !do iAt = 1, nAt
    !   write(18,'("Atom #: ",i3,"   position: ",3(f7.3,1x))') iAt,aBas(:,iAt)
    !enddo

    !write(*,'("Diagonal entries of HNF: ",3(i3,1x))') diag
    !write(*,'("Left transform:",/,3(3i3,1x,/),/)') transpose(L)
    !write(*,'("parent lattice vectors (columns):",/,3(3(f7.3,1x),/),/)') transpose(A)
    !write(*,'("parLatt inverse (columns):",/,3(3(f7.3,1x),/),/)') transpose(Ainv)
    !do iAt = 1, nAt
    !   write(*,'("Atom #: ",i3,"   position: ",3(f7.3,1x))') iAt,aBas(:,iAt)
    !enddo

    g = 0
    do iAt = 1, nAt ! Map each real space vector (atom position) into
       ! the group. We need to account for the fact that some of the
       ! atomic positions are offset by one of the d-vectors. So loop
       ! over the d-vectors and make sure that there is one mapping
       ! into the group that has no (negligible) non-integer parts. If
       ! not, something is wrong.
       mapped = .false.
       do iD = 1, nD
          gtemp = matmul(matmul(L,Ainv),aBas(:,iAt)-dset(:,iD))
          !write(*,'("shifted site vector",3(f8.4,1x))') aBas(:,iAt)-dset(:,iD)
          !write(*,'("iAt: ",i2,"    iD: ",i2,5x,3(f8.4,1x))') iAt, iD, gtemp
          if (equal(gtemp,nint(gtemp),eps)) then;
             mapped = .true.; dmember(iAt) = iD; exit
          endif
       enddo
       if (.not. mapped) stop "One of the atomic sites in find_labeling_from_atom_basis didn't map into the group"
       g(:,iAt) = nint(gtemp)
       g(:,iAt) = modulo(g(:,iAt),diag) ! Move each vector into first cell in g-space
    enddo

    ! Print out g-space representation of each vector
    !write(18,'("group list:")')
    !do i = 1, 3
    !   write(18,'(5x,200(i2,1x))') g(i,:)
    !enddo
    !write(18,'(5x,200(i2,1x))') dmember
    
    ! Find the original (unpermuted) group
    call make_member_list(diag,p)
    !write(18,'("original member list:")')
    !do i = 1, 3
    !   write(18,'(200(i2,1x))') p(i,:)
    !enddo
    allocate(p_with_d(4,nAt))
    
    ! Load up a 4xn*nD list of the group members
    n = nAt/nD ! volume factor of the supercell
    do iD = 1, nD ! Loop over each d-vector
       pLindx = (iD-1)*n+1 ! Index with a stride equal to volume factor
       ! need this because we want to set the group members for each d-member
       p_with_d(1:3,pLindx:pLindx+n-1) = p
       p_with_d(4,  pLindx:pLindx+n-1) = iD
    enddo
    !do i = 1, 4
    !   write(17,'(400(i2,1x))') p_with_d(i,:)
    !enddo
    
    call matrix_inverse(matmul(L,Ainv),T)
    !write(18,'("Map from group:",/,3(3(f7.3,1x),/),/)') transpose(T)
    
    ! If we want to do this, we need the supercell vectors
    !do iAt = 1, nAt
    !   postemp = matmul(T,p_with_d(1:3,iAt))+dset(:,p_with_d(4,iAt))
    !   call bring_into_cell(postemp
    !   write(18,'("Atom ",i3,"   from group (real space coords) ",3(f7.3,1x))') &
    !        iAt,matmul(T,p_with_d(1:3,iAt))+dset(:,p_with_d(4,iAt))
    !enddo
    
    ! Find the permutation of the original group that gives the input atom types
    
    call find_permutation_of_group_and_dset(p_with_d,g,dmember,perm)
    allocate(labeling(nAt))
    labeling = aTyp(perm)
    !write(18,'("Input atom labels: ",200(i2,1x))') aTyp
    !write(18,'("Group order:       ",200(i2,1x))') perm
    !write(18,'("Labeling:          ",200(i2,1x))') labeling
    
  ENDSUBROUTINE find_labeling_from_atom_basis
  
  !!<summary>Takes a 4xn list of the original group members and maps
  !!them onto the group representation of the current atomic
  !!basis---this generates the permutation of the original group that
  !!represents the current structure.</summary>
  !!<parameter name="g" regular="true">Unpermuted group.</parameter>
  !!<parameter name="gp" regular="true">Permuted group.</parameter>
  !!<parameter name="d" regular="true"></parameter>
  !!<parameter name="perm" regular="true">Permutation of gp</parameter>
  SUBROUTINE find_permutation_of_group_and_dset(g,gp,d,perm)
    integer, intent(in), dimension(:,:) :: g, gp ! unpermuted and permuted groups
    integer, intent(in), dimension(:)   :: d
    integer, intent(out) :: perm(:) ! permutation of gp

    integer n ! number of elements in the group (index of the superlattice * d-set size)
    logical skip(size(g,2))
    integer im, jm
    n = size(g,2); perm = 0
    skip = .false. ! This is just for efficiency
    do im = 1, n
       do jm = 1, n
          if (skip(jm)) cycle ! This is just for efficiency
          if (all(gp(:,jm)==g(1:3,im)).and.g(4,im)==d(jm)) then
             perm(im) = jm
             skip(jm) = .true.
             exit ! don't keep looking if you already found the match
          endif
       enddo ! jm
    enddo ! im
    !write (*,'(30i3)') perm
    if (any(perm==0)) then
       write(*,'(5x,"perm: ",200(i2,1x))') perm
       write(*,'(A)') "Maybe your epsilon is too large, and so some d-vectors appear to be the same after mapping in find_labeling_from_atom_basis"
       stop "mapping failed in find_permutation_of_group_and_dset"
    endif
  ENDSUBROUTINE find_permutation_of_group_and_dset

  !!<summary>This routine takes a derivative structure and finds its
  !!g-space representation. The routine returns all of the equivalent
  !!labelings of the superlattice so that these can be used to compare
  !!against possible matches in struct_enum.out.</summary>
  !!<parameter name="pLV" regular="true">Parent lattice vectors.</parameter>
  !!<parameter name="dset">"atomic basis" vectors of the
  !!multilattice.</parameter>
  !!<parameter name="sLV" regular="true">Superlattice vectors</parameter>
  !!<parameter name="aBas">Atom basis vectors</parameter>
  !!<parameter name="aTyp" regular="true">Labels for atoms in the
  !!basis of the superlattice.</parameter>
  !!<parameter name="HNF" regular="true">HNFs to generate the
  !!superlattice from pLV</parameter>
  !!<parameter name="Latdim" regular="true"></parameter>
  !!<parameter name="pLabel">The list of permuted labels (all
  !!equivalent).</parameter>
  !!<parameter name="eps" regular="true">Finite precision tolerance.</parameter>
  SUBROUTINE get_gspace_representation(pLV,dset,sLV,aBas,aTyp,HNF,LatDim,pLabel,eps)
    real(dp), intent(in) :: pLV(3,3) 
    real(dp), pointer ::  dset(:,:)  ! (in) 
    real(dp), intent(in) :: sLV(3,3) 
    real(dp), pointer :: aBas(:,:)   
    integer, intent(in) :: aTyp(:)   
    integer, intent(in) :: HNF(:,:,:)
    real(dp), intent(in) :: eps
    integer, intent(in) :: LatDim
    integer, pointer :: pLabel(:,:)  
    integer, pointer :: rotProdLabeling(:) ! List of which rotations
    ! producing each labeling permutation
    !logical err 

    integer, pointer :: aTypTemp(:), SNFlabel(:)
    real(dp), pointer :: aBasTemp(:,:), sLVlist(:,:,:)
    real(dp), dimension(3,3) :: pLVtemp, sLVtemp, parLattTest, T
    integer nAt, nuq, nP, SNF(3,3), iAt, ip, j, iop, iuq
    integer,pointer,dimension(:,:,:) :: HNFin, HNFout, L, R, SNFlist, uqSNF
    type(RotPermList) :: dRotList ! This is a list of permutations for the d-set 
    ! (needed as in put for several routines)
    type(opList), pointer :: fixOp(:) ! List of symops that fix the superlattice
    type(RotPermList), pointer :: LattRotList(:) ! List of rotation perms for the superlattice
    integer, pointer :: labeling(:), hnf_degen(:)

    integer, pointer :: label(:,:), digit(:)
    !debug
    !real(dp) Ainv(3,3)
    logical err
    !integer diag(3)
    ! debug
    
    !write(*,'("d-set: ",/,200(3(f7.3,1x),/))') (dset(:,iAt),iAt=1,size(dset,2))
    allocate(label(1,size(dset,2)))
    allocate(digit(size(dset,2)))
    label = 1
    digit = 1
    
    open(17,file="debug_gspace_rep.out")
    nAt = size(aBas,2)
    write(17,'("Number of atoms: ",i3)') nAt
    if(nAt/=size(aTyp)) stop "Input to get_gspace_representation is inconsistent: nAt"
    
    
    allocate(HNFin(3,3,1))
    HNFin(:,:,1) = HNF(:,:,1) ! The first HNF in the list is the one that is unrotated
    
    
    allocate(aTypTemp(nAt),aBasTemp(3,nAt))
    aTypTemp = aTyp; aBasTemp = aBas
    sLVtemp = sLV 
    
    call make_primitive(sLVtemp,aTypTemp,aBasTemp,.false.,eps)
    if(nAt/=size(aTypTemp)) then;
       write(17,'(/,"ERROR: The input structure wasn''t primitive")')
       write(17,'("atom type: ",20(i2,1x))') aTypTemp(:)
       write(17,'("number of atoms: ",i2,5x," size of aTyp:",i2)') nAt, size(aTypTemp)
       !stop "ERROR: input structure for get_gspace_representation was not primitive";endif
       write(*,'(A)') "WARNING: input structure for get_gspace_representation was not primitive"
    endif
    aTypTemp = 1; pLVtemp = sLV
    call make_primitive(pLVtemp,aTypTemp,aBasTemp,.false.,eps)
    write(17,'("Parent lattice (columns): ",/,3(3(f7.3,1x),/))') (pLVtemp(iAt,:),iAt=1,3)
    write(17,'("Parent lattice (in): ",/,3(3(f7.3,1x),/))') (pLV(iAt,:),iAt=1,3)
    write(17,'("Superlattice (columns): ",/,3(3(f7.3,1x),/))') (sLV(iAt,:),iAt=1,3)
    call matrix_inverse(pLVtemp,parLattTest,err)
    if(err) stop "Parent lattice of input superstructure is wrong"
    T = matmul(parLattTest,pLV)
    if(.not. equal(T,nint(T),eps)) stop "Input for get_gspace_representation is inconsistent"
    write(17,'("Size of supercell: ",i3)') abs(nint(determinant(sLV)/determinant(pLV)))
    write(17,'("d-set: ",/,200(3(f7.3,1x),/))') (dset(:,iAt),iAt=1,size(dset,2))
    !write(*,'("d-set: ",/,200(3(f7.3,1x),/))') (dset(:,iAt),iAt=1,size(dset,2))

    !** Calls to enumlib routines **

    ! This call generates a list of permutations for the d-set under symmetry operations
    ! of the parent lattice (need this for d-g table permutations)

    ! Wiley ! If gfortran compile breaks it's because we uncommented
    ! the write statement for the dRotList%RotIndx(:). Line 1071
    call get_dvector_permutations(pLV,dset,dRotList,LatDim,eps)
    write(17,'("d-set permutations: ",200(i3,1x))') dRotList%RotIndx(:)
    ! This call returns a list of operations that fix the
    ! superlattice. The routine expects a *list* of HNF matrices, but
    ! here we only need to pass in one because every one in the list is
    ! rotationally-equivalent. So the input and output lists are only
    ! one element long (in the last index). E.g., HNFin is a 3x3x1 array
    ! of integers (original lattice HNF)
    call remove_duplicate_lattices(HNFin,LatDim,pLV,dset,dRotList,HNFout,fixOp,LattRotList,sLVlist,hnf_degen,eps)
    write(17,'("Number of symmetry operations that fix the superlattice: ",i3,/)') size(fixOp(1)%rot,3)
    write(17,'(200(3(3f7.3,1x,/),"shift:",3(f7.3,1x),//))') &
         ((fixOp(1)%rot(j,:,iOp),j=1,3),fixOp(1)%shift(:,iOp),iOp=1,size(fixOp(1)%rot,3))
    ! This routine gets the SNF form of the HNF, returns a list of permutations effected by the
    ! rotation symmetries,
    call get_SNF(HNFout,L,SNFlist,R,LattRotList,uqSNF,SNFlabel,fixOp)
    ! Debug
    !print *,"HNF check",all(HNFin(:,:,1)==HNFout(:,:,1))
    !call matrix_inverse(pLV,Ainv,err)
    !diag = (/SNFlist(1,1,1),SNFlist(2,2,1),SNFlist(3,3,1)/)
    !write(*,'("Left transform:",/,3(3i2,1x,/))') transpose(L(:,:,1))
    !do iAt = 1,nAt
    !
    !   write(*,'("from group: # ",i3,5x," pos.:",3(i2,1x))') iAt, modulo(nint(matmul(L(:,:,1)&
    !        &,matmul(Ainv,aBas(:,iAt)))),diag)
    !enddo

    write(17,'("Left transform:",/,3(3i2,1x,/))') transpose(L(:,:,1))
    write(17,'("Smith Normal form:",/,3(3i2,1x,/))') transpose(SNFlist(:,:,1))
    write(17,'("Permutations:",/,8(24(i3,1x),/))') LattRotList(1)%RotIndx(:)
    
    SNF = SNFlist(:,:,1)
    
    ! This loads up the "perm" element of LattRotList
    call get_rotation_perms_lists(pLV,HNFout,L,SNFlist,fixOp,LattRotList,dRotList,eps)
    write(17,'("Rots Indx:",/,8(24(i3,1x),/))') LattRotList(1)%RotIndx(:)
    write(17,'("Permutation group (trans+rot):")')
    nP = size(LattRotList(1)%perm,1)
    do ip = 1, nP
       write(17,'("Perm #",i3,":",1x,200(i2,1x))') ip,LattRotList(1)%perm(ip,:)
    enddo
    
    ! debug
    !print *,"HNF check",all(HNFin(:,:,1)==HNFout(:,:,1))
    !call matrix_inverse(pLV,Ainv,err)
    !diag = (/SNFlist(1,1,1),SNFlist(2,2,1),SNFlist(3,3,1)/)
    !write(*,'("Left transform:",/,3(3i2,1x,/))') transpose(L(:,:,1))
    !do iAt = 1,nAt
    !
    !   write(*,'("from group: # ",i3,5x," pos.:",3(i2,1x))') iAt, modulo(nint(matmul(L(:,:,1)&
    !        &,matmul(Ainv,aBas(:,iAt)))),diag)
    !enddo


    call find_labeling_from_atom_basis(L(:,:,1),pLV,aBas,dset,aTyp,SNF,eps,labeling)
    write(17,'("Input atom labels: ",200(i2,1x))') aTyp
    write(17,'("Labeling:          ",200(i2,1x))') labeling
    
    ! Now that we know the labeling for the input atomic basis, find
    ! all possible permutations of that labeling. These can be
    ! compared to the entries in struct_enum.out
    ! Use the permutations effected by the rotations that fix the
    ! superlattice to generate labelings that are equivalent. The list
    ! of equivalent labelings will be used when we look for a match in
    ! the struct_enum file
    call find_equivalent_labelings(labeling,LattRotList,pLabel,rotProdLabeling)
    
    nuq = size(pLabel,1)
    write(17,'(/,"Number of unique labelings: ",i5)') nuq
    do iuq = 1, nuq
       write(17,'("Rot #:",i4,3x,"uq Labeling # :",i3,5x,"labeling:",1x,200(i1,1x))')&
            & rotProdLabeling(iuq),iuq,pLabel(iuq,:)
    enddo
    close(17)

  END SUBROUTINE get_gspace_representation

  !!<summary></summary>
  !!<parameter name="labeling" regular="true"></parameter>
  !!<parameter name="rotPerm" regular="true"></parameter>
  !!<parameter name="lab"></parameter>
  !!<parameter name="rotProdLab">Keep a list of which rotation created
  !!each labeling permutation</parameter>
  SUBROUTINE find_equivalent_labelings(labeling,rotPerm,lab,rotProdLab)
    integer, intent(in)           :: labeling(:)
    type(RotPermList), intent(in) :: rotPerm(:)
    integer, pointer              :: lab(:,:)
    integer, pointer              :: rotProdLab(:) 
    
    integer iuq, nuq, ip, nP, nAt
    integer, pointer              :: tempLabeling(:,:), tempRotProdLab(:)
    logical unique
    nP = size(rotPerm(1)%perm,1)
    nAt = size(labeling)

    iuq = 0; nuq = 0; allocate(tempLabeling(nP,nAt)); tempLabeling = 0
    allocate(tempRotProdLab(nP)); tempRotProdLab = 0
    do ip = 1, nP
       unique = .true.
       do iuq = 1, nuq
          if(all(labeling(rotPerm(1)%perm(ip,:))==tempLabeling(iuq,:))) then
             unique = .false.
             tempRotProdLab(iuq) = iP
             exit
          endif
       enddo
       if (unique) then
          nuq = nuq + 1
          tempLabeling(nuq,:) = labeling(rotPerm(1)%perm(ip,:))
       endif
    enddo
    allocate(lab(nuq,nAt))
    allocate(rotProdLab(nuq))
    lab = tempLabeling(1:nuq,:)
    rotProdLab = tempRotProdLab(1:nuq)
  END SUBROUTINE find_equivalent_labelings
  
  !!<summary>Reads in structure info from a struct_enum.out-type file
  !!and compares to the g-space representation of a test
  !!structure.</summary>
  !!<parameter name="fname" regular="true"></parameter>
  !!<parameter name="pLV" regular="true"></parameter>
  !!<parameter name="dset"></parameter>
  !!<parameter name="HNFin">List of HNFs (all equivalent).</parameter>
  !!<parameter name="SNF" regular="true"></parameter>
  !!<parameter name="LatDim" regular="true"></parameter>
  !!<parameter name="pLabel">Permuted labelings</parameter>
  !!<parameter name="match" regular="true"></parameter>
  !!<parameter name="eps" regular="true">Finite precision tolerance.</parameter>
  SUBROUTINE find_match_in_structenumout(fname,pLV,dset,HNFin,SNF,LatDim,pLabel,match,eps) 
    character(80), intent(in) :: fname
    real(dp), intent(in), dimension(3,3) :: pLV
    integer, intent(in), dimension(3,3) :: SNF
    real(dp), pointer :: dset(:,:) ! (in)
    integer, pointer :: pLabel(:,:) ! (in)  
    integer, pointer :: HNFin(:,:,:) ! (in) 
    integer, intent(in) :: LatDim
    integer, intent(out) :: match
    real(dp) :: eps

    character(80) title, bulksurf, dummy
    character(maxLabLength) :: labeling
    real(dp) :: p(3,3), Ainv(3,3)
    real(dp), allocatable :: dvec(:,:)
    integer, allocatable :: ilabeling(:)
    integer iStr, volume, HNF(3,3), iLab, nLab, ioerr, i, nD, k, strN, hnfN, sizeN
    integer a, b, c, d, e, f, pgOps, diag(3), L(3,3), n, iHNF, nHNF
    logical foundLab, foundHNF
    
    open(13,file="readcheck_find_match_in.out")
    write(13,'(5(/),"This file echoes the input for the find_match_in_structenumout routine and")')
    write(13,'("prints some of the computed quantities as well. It is intended to")')
    write(13,'("be useful for debugging and fixing bad input")')
    write(13,'("***************************************************")')
    
    volume = determinant(SNF)
    !call read_nth_line_from_enumlist()
    open(11,file=fname,status='old',iostat=ioerr)
    if(ioerr/=0)then; write(*,'("Input file doesn''t exist:",a80)') trim(fname);stop;endif
    ! Read in the title from the struct_enum.out file
    read(11,'(a80)') title; title = trim(title)
    write(13,'("Title: ",a80)') title
       
    ! Read in surf/bulk mode marker
    read(11,'(a1)') bulksurf
    write(13,'("bulk/surf: ",a1)') bulksurf
    call ucase(bulksurf)
    if((LatDim==2 .and. bulksurf=="B") .or. (LatDim==3 .and. bulksurf=="S")) then 
       print *,"The bulk/surf parameter in "//trim(fname)//" not consistent"
       write(*,'("Bulk/surf parameter is: ",a1)') bulksurf
       write(*,'("LatDim input in call is: ",i1)') LatDim
       stop
    endif
    
    ! Read in the parent lattice vectors
    do i = 1,3; read(11,*) p(:,i); enddo
    call matrix_inverse(p,Ainv)
    !write(13,'("Parent lattice:",/,3(3f7.3,1x,/))') (p(:,i),i=1,3) 
    !write(13,'("Parent lattice inverse matrix:",/,3(3f7.3,1x,/))') (Ainv(:,i),i=1,3) 
    if (.not. equal(determinant(pLV),determinant(p),eps)) then
       write(13,'("Volume of input structure parent: ",f7.3)') determinant(pLV)
       write(13,'("Volume of struct_enum parent: ",f7.3)') determinant(p)
       !stop "Volumes of parent lattices don't match"
    endif
    !if (.not. equal(matmul(Ainv,pLV),nint(matmul(Ainv,pLV)),eps)) stop "Parent lattices don't match"
    
    ! Read in the number of d-vectors, then read each one in
    read(11,*) nD
    write(13,'("Number of d-vectors: ",i3)') nD
    if(nD/=size(dset,2)) stop "Numbers of d-vectors don't match"
    allocate(dvec(3,nD),ilabeling(nD*volume))
    do i = 1,nD; read(11,*) dvec(:,i); enddo
    write(13,'("d-vectors:",/,40(3f8.4,1x,/))') (dvec(:,i),i=1,nD) 
       

    ! Read in the number of labels, i.e., binary, ternary, etc.
    read(11,'(i2)') k
    write(13,'("Number of labels (k): ",i2)') k 
    
    read(11,*) ! Skip starting and ending cell sizes
    read(11,*) eps
    write(13,'("Finite precision parameter (epsilon): ",g17.10)') eps 
    
    do 
       read(11,*) dummy
       if(dummy=="start") exit
    enddo
    write(13,'("Found ""start"" label. Now reading in structures line-by-line ")')  

    nHNF = size(HNFin,3)
    
    ! Read in the info for each structure in the file
    write(13,'("Looking for a structure of size (volume factor): ",i3)') volume
    write(13,'("strN, hnfN, sizeN, n, pgOps, diag, a,b,c,d,e,f, L, labeling")')  
    iStr = 0; match = 0
    nLab = size(pLabel,1)
    do 
       iStr = iStr + 1
       read(11,*,iostat=ioerr) strN, hnfN, sizeN, n, pgOps, diag, a,b,c,d,e,f, L, labeling
       if(ioerr/=0) exit
       read(labeling,'(500i1)') ilabeling(1:n*nD)
       L = transpose(L) ! Listed with columns as the fast index in the
       !struct_enum.out but reads in with rows as the fast index.
       write(13,'(i11,1x,i7,1x,i8,1x,i2,1x,i2,1x,3(i3,1x),1x,6(i3,1x),1x,1x,50i1)') &
            strN, hnfN, sizeN, n, pgOps, diag, a,b,c,d,e,f,ilabeling
       ! Define the full HNF matrix
       HNF = 0
       HNF(1,1) = a; HNF(2,1) = b; HNF(2,2) = c
       HNF(3,1) = d; HNF(3,2) = e; HNF(3,3) = f
       if(n/=volume) then
          write(13,'("Volume doesn''t match for str #:",i8," ---Skipping")')  strN
          cycle
       endif
       foundHNF = .false.
       do iHNF = 1, nHNF
          if(.not. all(HNF==HNFin(:,:,iHNF))) cycle
          foundHNF = .true.
       enddo
       if(.not. foundHNF) then
          write(13,'("HNF doesn''t match for str #:",i8," ---Skipping")') StrN
          cycle
       endif
       foundLab = .false.
       do iLab = 1, nLab
          if(all(ilabeling==pLabel(iLab,:))) then
             foundLab = .true.
             if (match/=0) then
                close(13)
                print *, "The first two matching structures are:"
                print *,match,iStr
                stop "BUG! Found more than one match in struct_enum.out file"
             endif
             match = iStr
             exit
          endif
       enddo
       if(.not. foundLab) then
          write(13,'("Labeling didn''t match for str #:",i8)') strN
       else
          write(13,'("Structure number:",i8," is a match!")') strN
       endif
    enddo
    close(13)
    deallocate(ilabeling)
    
  End SUBROUTINE find_match_in_structenumout

  !!<summary>This routine is based on find_match_in_structenumout but
  !!relies on the supply of two structures only.</summary>
  !!<parameter name="LatDim" regular="true"></parameter>
  !!<parameter name="pLV" regular="true"></parameter>
  !!<parameter name="dset"></parameter>
  !!<parameter name="eps" regular="true">Finite precision tolerance.</parameter>
  !!<parameter name="A_HNF" regular="true">HNF of stucture A.</parameter>
  !!<parameter name="A_labeling" regular="true">Labeling of structure A.</parameter>
  !!<parameter name="B_HNFlist" regular="true">{ entry, entry, HNF# }
  !!List of HNFs (all equivalent)</parameter>
  !!<parameter name="B_labelingList" regular="true">{ labeling#, entry
  !!} permuted labelings</parameter>
  !!<parameter name="match" regular="true">Returns true if structures match.</parameter>
  SUBROUTINE compare_two_gstructures(LatDim,pLV,dset,eps, A_HNF, A_labeling, B_HNFlist, B_labelingList, match) 
    real(dp), intent(in), dimension(3,3) :: pLV
    integer, intent(in)                  :: LatDim
    real(dp), pointer                    :: dset(:,:)     ! (in)
    ! "original" structure A:
    integer, intent(in) :: A_HNF(:,:)          ! (in) HNF
    integer, intent(in) :: A_labeling(:)       ! (in) labeling
    ! structure B that is to be compared to the "original" structure:
    integer, intent(in) :: B_labelingList(:,:) ! (in) 
    integer, intent(in) :: B_HNFlist(:,:,:)    ! (in) 
    ! do both structures match?
    logical, intent(out) :: match
    
    real(dp) :: eps

    real(dp) :: p(3,3), Ainv(3,3)
    integer iLab, nLab
    integer iHNF, nHNF
    logical foundLab, foundHNF
    
    call matrix_inverse(p,Ainv)
    
    nHNF = size(B_HNFlist,3)
    nLab = size(B_labelingList,1)
    
    match = .false.
    
    foundHNF = .false.
    do iHNF = 1, nHNF
       if(.not. all(A_HNF==B_HNFlist(:,:,iHNF))) cycle
       foundHNF = .true.
       exit
    enddo
    if(.not. foundHNF) then  
       return
    endif
    foundLab = .false.
    do iLab = 1, nLab
       if(.not. all(A_labeling==B_labelingList(iLab,:))) cycle
       foundLab = .true.
       exit
    enddo
    
    if (foundHNF .and. foundLab) match = .true.
    return
  End SUBROUTINE compare_two_gstructures

END MODULE enumeration_utilities
