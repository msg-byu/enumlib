!!<summary>Contains the subroutines needed for the addition of the
!!arrows to the tree class enumeration.</summary>
Module arrow_related
use enumeration_types

implicit none
private
public  generateIndexFromArrowing, generateArrowingFromIndex, arrow_concs, write_arrow_labeling
contains
  
  !!<summary>Finds the integer corresponding to the arrow array.</summary>
  !!<parameter name="arrowing" regular="true">An integer array of the arrow
  !!directions.</parameter>
  !!<parameter name="dim" regular="true">The number of directions the arrows can
  !!point.</parameter>
  !!<parameter name="arrow_index" regular="true">The output integer index.</parameter>
  subroutine generateIndexFromArrowing(arrowing, dim, arrow_index)
    integer, intent(in) :: arrowing(:)
    integer, intent(in) :: dim
    integer, intent(out) :: arrow_index

    !!<local name="temp_arrowing">A temporary copy of the arrowing.</local>
    integer :: temp_arrowing(size(arrowing))
    integer :: i
    
    temp_arrowing = arrowing - 1

    arrow_index = 0
    do i=1, size(arrowing)
       arrow_index = arrow_index + temp_arrowing(i)*dim**(i-1)
    end do    
  end subroutine generateIndexFromArrowing
  
  !!<summary>Turns an index into an arrowing.</summary>
  !!<parameter name="index" regular="true">The integer index for the arrowing.</parameter>
  !!<parameter name="nArrows" regular="true">The number of arrows.</parameter>
  !!<parameter name="dim" regular="true">The number of directions the arrows can
  !!point.</parameter>
  !!<parameter name="arrowing" regular="true">The integer array of the arrow
  !!directions.</parameter>
  subroutine generateArrowingFromIndex(index, nArrows, dim, arrowing)
    integer, intent(in) :: index, nArrows, dim
    integer, intent(out) :: arrowing(nArrows)

    !!<local name="base">The arrow divisor.</local>
    !!<local name="temp_index">A copy of the index.</local>
    integer :: base, temp_index
    integer :: i
    
    temp_index = index
    arrowing = 0
    do i = 1, nArrows
       base = dim**(nArrows-i)
       arrowing(nArrows+1-i) = temp_index/base
       temp_index = temp_index - base*arrowing(nArrows+1-i)
    end do
    arrowing = arrowing + 1
  end subroutine generateArrowingFromIndex

  !!<summary>Finds the new concentration for the system with arrows
  !!included as well as the mapping that will restore the original
  !!atomic species information when a unique arrangement has been
  !!identified."></summary>
  !!<parameter name="conc" regular="true">The integer array of concentrations.</parameter>
  !!<parameter name="arrows" regular="true">The number of arrows for each atomic
  !!species.</parameter>
  !!<parameter name="a_conc" regular="true">The resultant concentrations with the arrows
  !!included.</parameter>
  !!<parameter name="conc_map" regular="true">A mapping from the new concentrations
  !!to the old.</parameter>
  subroutine arrow_concs(conc,arrows,a_conc,conc_map)
    integer, intent(in) :: conc(:)
    integer, intent(in) :: arrows(:)
    integer, allocatable, intent(out) :: a_conc(:)
    integer, allocatable, intent(out) :: conc_map(:,:)

    !!<local name="new_species">A temporary array of the new species
    !!to be included in the concentration list.</local>
    !!<local name="new_concs">An array of the concentrations of atoms
    !!not being displaced.</local>
    !!<local name="temp_map">A temporary array for the concentration map.</local>
    !!<local name="new_size">The size of the new concentration array.</local>
    integer :: new_species(size(conc)), new_concs(size(conc))
    integer :: temp_map(size(conc),2)
    integer :: i, j, new_size

    temp_map = 0
    ! First we need to find out the concentrations of the arrowed and
    ! unarrowed atoms of each species
    j = 1
    do i=1, size(conc)
       if ((conc(i) > 0) .and. (arrows(i) > 0)) then
          if (arrows(i) <= conc(i)) then
             new_concs(i) = conc(i) - arrows(i)
             new_species(j) = arrows(i)
          else
             new_concs(i) = 0
             new_species(j) = conc(i)
          end if
          temp_map(j,1) = size(conc) + j
          temp_map(j,2) = i
          j = j + 1
       else
          new_concs(i) = conc(i)
       end if
    end do
    j = j - 1

    new_size = size(conc) + j
    allocate(a_conc(new_size),conc_map(j,2))
    ! Now fill in the a_conc array
    do i = 1, new_size
       if (i <= size(conc)) then
          a_conc(i) = new_concs(i)
       else
          a_conc(i) = new_species(i-size(conc))
       end if
    end do

    !Now fill in the complete concentration mapping
    do i=1, j
       conc_map(i,:) = temp_map(i,:)
    end do
  end subroutine arrow_concs
  
  !!<summary>Writes a single labeling to file for the enum4
  !!code.</summary>
  !!<parameter name="n" regular="true">The size of the supercell.</parameter>
  !!<parameter name="Tcnt" regular="true">Counter for the total number
  !!of labelings.</parameter>
  !!<parameter name="Scnt" regular="true">Counter for the number of
  !!labelings for this size.</parameter>
  !!<parameter name="Hcnt" regular="true">Counter for the HNFs.</parameter>
  !!<parameter name="HNFi" regular="true">Index in the permIndx
  !!corresponding to the current block of the HNFs.</parameter>
  !!<parameter name="hnf_degen" regular="true">The degeneracy of the HNFs.</parameter>
  !!<parameter name="fixOp" regular="true">Lattice fixing operations (type opList).</parameter>
  !!<parameter name="SNFlist" regular="true">List of the SNFs.</parameter>
  !!<parameter name="HNFlist" regular="true">List of the HNFs.</parameter>
  !!<parameter name="L" regular="true">List of the left transforms.</parameter>
  !!<parameter name="permIndx" regular="true">List of the different permutations
  !!groups.</parameter>
  !!<parameter name="equivalencies" regular="true">The list of
  !!equivalencies of the system.</parameter>
  !!<parameter name="labeling" regular="true">The labeling to be written to file.</parameter>
  !!<parameter name="arrow_label" regular="true">The arrows to b written to file.</parameter>
  subroutine write_arrow_labeling(labeling,arrow_label,n,Tcnt,Scnt,Hcnt,HNFi,hnf_degen,fixOp,SNFlist,HNFlist,L,equivalencies,permIndx)
    integer, intent(in)      :: n, Hcnt, HNFi
    integer, intent(inout)   :: Scnt, Tcnt
    integer, intent(in)      :: SNFlist(:,:,:), HNFlist(:,:,:), L(:,:,:)
    integer, intent(in)      :: equivalencies(:), hnf_degen(:), permIndx(:), labeling(:), arrow_label(:)
    type(opList), intent(in) :: fixOp(:)
    
    !!<local name="lab_degen">The degeneracy of this label.</local>
    !!<local name="iHNF">Counter that loops over the HNFs</local>
    !!<local name="vsH">Vector for matching the HNF to the list.</local>
    !!<local name="struct_enum_out_formatstring">Format for output file.</local>
    !!<local name="dummy">Dummy string.</local>
    !!<local name="status">Allocation exit flag</local>
    !!<local name="i">Loop variable</local>
    !!<local name="jHNF">Which HNF this is in the permIndx</local>
    !!<local name="nHNF">How many of this HNF are there.</local>
    integer :: lab_degen, iHNF, status, i, jHNF, nHNF
    integer, allocatable :: vsH(:)
    character(105) :: struct_enum_out_formatstring
    character(3) :: dummy

    ! Labeling Postprocessing data
    !!<local name="nALLD">Size of full dset.</local>
    !!<local name="allD">help array: (/1,2,3,....,nALLD/)</local>
    !!<local name="pplabeling">Labeling after postprocessing.</local>
    !!<local name="postprocessLabeling">True if we need to perform
    !!postprocessing.</local>
    !!<local name="allD2LabelD">Respective d-vector ID in the labeling dset.</local>
    integer              :: nAllD
    integer, allocatable :: allD(:)
    integer, allocatable :: pplabeling(:)
    logical :: postprocessLabeling
    integer, allocatable :: allD2LabelD(:)


    ! Postprocessing labelings: setup
    nAllD = size(equivalencies)
    allocate(allD(nAllD)); allD = (/ (i,i=1,nAllD) /)
    allocate(allD2LabelD(nAllD));
    ! 1) Check whether we have to postprocess the labeling before writing
    !    it out Postprocessing is needed if we do not want to enumerate
    !    all dset members of a primitive unit cell due to some
    !    equivalencies.  For example, in a 1x1 symmetric surface slab (
    !    (*) denotes an atom ):
    !
    !      ------------------------------------------------ surface
    !         (*)  topmost surface layer, dvector# 1
    !         (*)  dvector# 2
    !         (*)  dvector# 3
    !         (*)  bottommost surface layer, dvector# 4
    !      ------------------------------------------------ surface
    !
    !    the topmost and the bottommost atom should always have the same
    !    occupancy, as well as dvector 2 and dvector 3 should. You can
    !    therefore specify the following equivalency list:
    !
    !         equivalency of dvector# | 1 2 3 4
    !         ---------------------------------   
    !         equivalency list        | 1 2 2 1
    !
    !    which means that dvector# 1 and dvector# 4 have to have the same
    !    occupancy, they are equivalent by enumeration. The same is true
    !    for dvector# 2 and dvector# 3
    !
    !    In this example, the enumeration code should only find
    !    enumerations of dvector# 1 and dvector# 2.  The occupations of
    !    dvector# 3 and dvector# 4 are then constructed in a
    !    postprocessing step.  The postprocessing step takes the
    !    enumerated form (i.e. dvectors# 1 and 2, e.g. a labeling 0101 for
    !    two unit cells) and tranforms it into a form that is valid for
    !    ALL dvectors, e.g. labeling 01010101 for two unit cells).
    !
    postprocessLabeling = .not. (all(  abs( equivalencies-allD ) ==0))
    allocate(pplabeling(n*nAllD))
    ! this is ok whether we do postprocessing (nAllD>nD) or not
    ! (nAllD==nD) 2) if yes: - prepare the postprocessed labeling
    ! (pplabeling) - generate a map: full dset -> enum dset that tells
    ! me for a dset member of the full dset what dset member in the
    ! enumeration dset it corresponds to
    if (postprocessLabeling) then
       allD2LabelD = (/(count(pack(allD,allD<=equivalencies(i))==equivalencies),i=1,nAllD)/)
       ! this construction is best explained by an example: suppose we
       ! have a dset 1,2,3,4,5 and the equivalent list is 1,4,4,4,5
       ! (so, dset member 1,4,5 will be enumerated, having positions
       ! 1,2,3 in labelings) the map should accomplish: 1 -> 1, 2 ->
       ! 2, 3 -> 2, 4 -> 2, 5 -> 3 first, the pack operation take only
       ! those in the dset <= equivalent point, then, the number of
       ! truly enumerated points is counted (truly enumerated points
       ! are points for which dset=equivalencies
    endif

    lab_degen = 0
    nHNF = count(permIndx==HNFi)
    allocate(vsH(nHNF),STAT=status); if (status/=0) stop "Allocation failed in write_single_labelings: vsH"
    
    ! Packing...
    vsH = pack((/(i,i=1,size(HNFlist,3))/), HNFi==permIndx); 

    write(dummy,'(I3)') n*nAllD
    
    struct_enum_out_formatstring = '(i11,1x,i9,1x,i7,1x,i8,1x,i8,1x,i11,1x,i3,2x,i4,2x,3(i2,1x),2x,6(i2,1x),2x,9(i4,1x),2x,'//trim(dummy)//'i1,4x,'//trim(dummy)//'i1)'

    if (postprocessLabeling) then
       ! see the comments at the beginning of the current routine
       call postprocess_labeling(n,nAllD,labeling,pplabeling,allD2LabelD) 
    else
       pplabeling = labeling ! nothing changes
    endif

    do iHNF = 1, nHNF ! Write this labeling for each corresponding HNF
       jHNF = vsH(iHNF) ! Index of matching HNFs
       Tcnt = Tcnt + 1; Scnt = Scnt + 1
       write(14,struct_enum_out_formatstring) &
            Tcnt, Hcnt+iHNF,hnf_degen(jHNF),lab_degen,lab_degen*hnf_degen(jHNF),&
            Scnt,n,size(fixOp(jHNF)%rot,3),SNFlist(1,1,jHNF),SNFlist(2,2,jHNF),&
            SNFlist(3,3,jHNF),HNFlist(1,1,jHNF),HNFlist(2,1,jHNF),HNFlist(2,2,jHNF),&
            HNFlist(3,1,jHNF),HNFlist(3,2,jHNF),HNFlist(3,3,jHNF),transpose(L(:,:,jHNF)),&
            pplabeling, arrow_label
    enddo ! loop over HNFs

  contains

    !!<summary>Purpose: take an enumeration labeling (for selected,
    !!non-equivalent (by enumeration) dvectors) and construct the full
    !!labeling for all dvectors. It makes sure that two dvectors in
    !!the same primitive unit cell of the parent lattice get the SAME
    !!labeling always. See also comments at the beginning of
    !!write_labelings.</summary>
    !!<parameter name="nUC" regular="true">number of unit
    !!cells.</parameter>
    !!<parameter name="nAllD" regular="true">number of d-vectors in
    !!the new labeling.</parameter>
    !!<parameter name="oldlabeling" regular="true">the old
    !!labeling</parameter>
    !!<parameter name="newlabeling" regular="true">the wanna-be new
    !!labeling.</parameter>
    !!<parameter name="newD2oldD" regular="true">{ new D# }: a map
    !!newD -> oldD</parameter>
    subroutine postprocess_labeling(nUC,nAllD,oldlabeling,newlabeling,newD2oldD)
      integer, intent(in) :: nUC, nAllD          
      integer, intent(in) :: oldlabeling(:)      
      integer, intent(out):: newlabeling(:)      
      integer, intent(in) :: newD2oldD(:)        

      !!<local name="newlab_pos">The position in the new label.</local>
      !!<local name="oldlab_pos">The position in the old label.</local>
      !!<local name="iD">Variable for looping.</local>
      !!<local name="iUC">Variable for looping.</local>
      integer :: newlab_pos, oldlab_pos
      integer :: iD,iUC
      
      do iD=1,nAllD
         do iUC=1,nUC
            newlab_pos = (iD-1)*nUC + iUC               ! position in the new labeling
            oldlab_pos = (newD2oldD(iD)-1)*nUC + iUC    ! corresponding
            ! position in the old labeling
            newlabeling(newlab_pos) = oldlabeling(oldlab_pos)
         enddo
      enddo

    end subroutine postprocess_labeling

  end subroutine write_arrow_labeling
end Module arrow_related
