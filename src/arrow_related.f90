!!<summary>Contains the subroutines needed for the addition of the
!!arrows to the tree class enumeration.</summary>
Module arrow_related
use enumeration_types

implicit none
private
public  generateIndexFromArrowing, generateArrowingFromIndex, arrow_concs
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
    integer :: arrow_i, arrow_j

    temp_arrowing = arrowing - 1

    arrow_index = 0
    arrow_j = 1
    do arrow_i=1, size(arrowing)
       if (temp_arrowing(arrow_i) >= 0) then
          arrow_index = arrow_index + temp_arrowing(arrow_i)*dim**(arrow_j-1)
          arrow_j = arrow_j + 1
       end if
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
    integer :: arrow_i
    
    temp_index = index
    arrowing = 0
    do arrow_i = 1, nArrows
       base = dim**(nArrows-arrow_i)
       arrowing(nArrows+1-arrow_i) = temp_index/base
       temp_index = temp_index - base*arrowing(nArrows+1-arrow_i)
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
  !!<parameter name="nArrows" regular="true">The number of sites with arrows.</parameter>
  subroutine arrow_concs(conc,arrows,a_conc,conc_map, nArrows)
    integer, intent(in) :: conc(:)
    integer, intent(in) :: arrows(:)
    integer, allocatable, intent(out) :: a_conc(:)
    integer, allocatable, intent(out) :: conc_map(:,:)
    integer, intent(out) :: nArrows

    !!<local name="new_species">A temporary array of the new species
    !!to be included in the concentration list.</local>
    !!<local name="new_concs">An array of the concentrations of atoms
    !!not being displaced.</local>
    !!<local name="temp_map">A temporary array for the concentration map.</local>
    !!<local name="new_size">The size of the new concentration array.</local>
    integer :: new_species(size(conc)), new_concs(size(conc))
    integer :: temp_map(size(conc),2)
    integer :: species_i, species_j, new_size

    nArrows = 0
    temp_map = 0
    ! First we need to find out the concentrations of the arrowed and
    ! unarrowed atoms of each species
    species_j = 1
    do species_i=1, size(conc)
       if ((conc(species_i) > 0) .and. (arrows(species_i) > 0)) then
          if (arrows(species_i) <= conc(species_i)) then
             new_concs(species_i) = conc(species_i) - arrows(species_i)
             new_species(species_j) = arrows(species_i)
             nArrows = nArrows + arrows(species_i)
          else
             new_concs(species_i) = 0
             new_species(species_j) = conc(species_i)
             nArrows = nArrows + conc(species_i)
          end if
          temp_map(species_j,1) = size(conc) + species_j
          temp_map(species_j,2) = species_i
          species_j = species_j + 1
       else
          new_concs(species_i) = conc(species_i)
       end if
    end do
    species_j = species_j - 1

    new_size = size(conc) + species_j
    allocate(a_conc(new_size),conc_map(species_j,2))
    ! Now fill in the a_conc array
    do species_i = 1, new_size
       if (species_i <= size(conc)) then
          a_conc(species_i) = new_concs(species_i)
       else
          a_conc(species_i) = new_species(species_i-size(conc))
       end if
    end do

    !Now fill in the complete concentration mapping
    conc_map(:,:) = temp_map(:species_j,:)
  end subroutine arrow_concs
  
end Module arrow_related
