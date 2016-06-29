!!<summary>AUTHORS: Conrad W. Rosenbrock, Wiley S. Morgan (October 2014)
!! Classes to support the calculation of coefficients for specific terms in a product
!! of multinomials. Construct a product class by specifying the exponent and target term
!! and then add multinomials using the Product instance's append(). The coefficient is
!! then available from the coeff().
!!</summary>
MODULE classes
use itertools
IMPLICIT NONE
public polya, lint
  !!<summary>Big integer -10^18..10^18 range</summary>
  integer, parameter:: lint=selected_int_kind(18)
  !!<summary>Represents an exponent-limited sequence with a single root. Here, sequence represents a
  !!sequence of integer values k_1, k_2...k_j that are the exponents of a single term in a multinomial.
  !!The root of the sequence is one of the k_i; its children become sets of sequences including
  !!variables to the right of i</summary>
  type, public :: Sequence
     !!<member name="root">The exponent of the variable to the left in the multinomial term.</member>
     !!<member name="parent">A Sequence instance for the variable to the *left* of this one 
     !!(i.e. has index i-1).</member>
     !!<member name="varcount">The number of variables in the multinomial (also length of targets).</member>
     !!<member name="kids"></member>
     !!<member name="used">The sum of exponents to the left of this variable sequence.</member>
     !!<member name="kidcount">The number of child sequences to the right of this sequence.</member>
     integer :: root, used, varcount, kidcount
     class(Sequence), pointer :: parent => null()
     class(Sequence), pointer :: kids(:)
     contains
       procedure, public :: initialize => sequence_initialize
       procedure, public :: expand
       procedure, public :: finalize => sequence_finalize
  end type Sequence

  !!<summary>Represents a multinomial expansion.</summary>
  type, public :: Multinomial 
     !!<member name="pawer">the powers on each of the *unexpanded* variables in the multinomial;
     !!of the form (x^2+y^2+z^2) => 2. No support is provided for coefficients on the
     !!variables inside the brackets (only Polya implementation).</member>
     !!<member name="exponent">the exponent of the entire multinomial.</member>
     !!<member name="powersum">The maximum allowed exponent on any term in the expansion.</member>
     !!<member name="possible_powers">The list of possible exponents a single variable can have in
     !!any of the terms in the expansion.</member>
     integer, allocatable :: power
     integer :: exponent, powersum
     integer, allocatable :: possible_powers(:)
     contains
       procedure, public :: initialize => multinomial_initialize
       procedure, public :: nchoosekm
  end type Multinomial

  !!<summary>Represents a product of multinomials for which only a single term is interesting.</summary>
  type, public :: Product
     !!<member name="coefficient">The scalar integer multiplying this product of multinomials.</member>
     !!<member name="targets">a list of exponents for the only interesting term in the product. The
     !!list is in the order that the variables appear in each multinomial.</member>
     !!<member name="multinoms">An array containing the multinomials.</member>   
     !!<member name="multinom_rd">A 2D array of the r-cycle and degeneracy values that define each
     !!multinomial.</member>
     integer :: coefficient
     integer, allocatable :: targets(:)
     class(Multinomial), allocatable :: multinoms(:)
     integer, allocatable :: multinom_rd(:,:)
     contains
       procedure, public :: initialize => product_initialize
       procedure, public :: coeff 
       procedure, public :: sum_sequence
  end type Product

contains
  !!<summary>Initializes a sequence collector for a variable. 'Term' refers to a product
  !!of variables like x^i.y^j.z^k, then x has index 0, y has 1, etc.</summary>
  !!<parameter name="root" regular="true">The exponent of the variable to the left 
  !!in the multinomial term.</parameter>
  !!<parameter name="possibles" regular="true">A list of possible values for each variable 
  !!in the multinomial.</parameter>
  !!<parameter name="depth" regular="true">The index of the variable being sequenced in 
  !!the term. Should be set to \emph{2} for the top level.</parameter>
  !!<parameter name="powersum" regular="true">The maximum value that the sum of exponents in the sequence 
  !!is allowed to have.</parameter>
  !!<parameter name="targets" regular="true">An array of the target exponents on each variable.</parameter>
  !!<parameter name="parent">A Sequence instance for the variable to the *left* of this one 
  !!(i.e. has index i-1).</parameter>
  recursive subroutine sequence_initialize(this, root, possibles, depth, powersum, targets, parent)
    class(Sequence) :: this
    integer, intent(in) :: powersum, root, depth
    integer, allocatable, intent(in) :: targets(:)
    integer, intent(in) :: possibles(:)
    class(Sequence), target, optional :: parent
    
    !!<local name="kidroots">An array of the exponents af the children, ie the possibles-root 
    !!that still meet the conditions</local>
    !!<local name="kidcount">The number of kids that each depth has to its right.</local>
    !!<local name="p">Temporary variable for analyzing each possible exponent for the
    !!variable at the current depth.</local>
    integer, allocatable :: kidroots(:)
    integer :: i, p, kidcount

    this%root = root
    if (present(parent)) then
       this%parent => parent
       this%used = this%root + this%parent%used
    else
       this%used = this%root
    end if

    this%varcount = size(targets, 1)
    allocate(kidroots(size(possibles, 1)))
    !We only keep recursively defining sequences until we run out of variables in
    !the term. Possibles is a list of possible exponents for each variable in the
    !term and has the same number of items as variables in the term.
    if (depth .le. this%varcount) then
       !Filter the possible values for the variable being considered based on the
       !exponent of the multinomial. When multinomials are expanded, the sum of
       !the exponents in any term must be less than the exponent on the multinomial
       !times the maximum power of any of its (unexpanded) terms.
       !We find all the possible values for this variable by ensuring that:
       !1) it's exponent is compatible with the exponents of all variables to the left of it.
       !2) the exponent we are suggesting is in the list of possible values for the variable.
       !3) the exponent remains positive.
       kidcount = 0
       do i = 1, size(possibles,1)
          p = possibles(i)
          if (p-root >= 0) then
             if ((p-root <= targets(depth)) .and. &
                  (abs(p-root) <= powersum-this%used) .and. (mod(abs(p-this%used), possibles(2)) .eq. 0)) then
                kidcount = kidcount+1
                kidroots(kidcount) = p-root
             end if
          end if
       end do
       allocate(this%kids(kidcount))
 
       !Construct possible exponents for the variable to the right of this one based on the
       !possible choices we just calculated. This gets done recursively until we have a tree
       !of possible exponent values for the variables from left to right.
       this%kidcount = 0
       do i = 1, kidcount
          call this%kids(i)%initialize(kidroots(i), possibles, depth+1, powersum, targets, this)
          this%kidcount = this%kidcount + this%kids(i)%kidcount
       end do

       if (this%kidcount .eq. 0) this%kidcount = kidcount
    else 
       this%kids => null()
       this%kidcount = 0
    end if
  end subroutine sequence_initialize
  
  !!<summary>Cleans up the local variables so that @CREF[this.initialize()] can be called on the
  !!*same* instance that is already allocated.</summary>
  subroutine sequence_finalize(this)
    class(Sequence) :: this
    deallocate(this%kids)
  end subroutine sequence_finalize

  !!<summary>Recursively generates a list of all relevant sequences for this multinomial term.</summary>
  !!<parameter name="sequences" regular="true">The result of the expansion; is allocated 
  !!within this routine. Returns an array with dimensions (sum(allkids), varcount).</parameter>
  !!<parameter name="pstart" regular="true">The row in 'sequences' to start filling in during the current 
  !!recursion. Should be 1 for the top level.</parameter>
  !!<parameter name="varindex" regular="true">The index of the variable in the polynomials. Should be 1
  !!for the top level.</parameter>
  recursive subroutine expand(this, sequences, pstart, varindex)
    class(Sequence) :: this 
    integer, intent(inout), allocatable :: sequences(:,:)
    integer, intent(in) :: pstart, varindex

    !!<local name="j,k">Iterators for the current sequences kids and the kids' variables
    !!respectively during this level of expansion.</local>
    !!<local name="cursor">Points to the current row in the 'sequences' result that is
    !!having its value updated.</local>
    !!<local name="start">Updatable copy of the 'pstart' parameter.</local>
    integer :: j, k, cursor
    integer :: start

    !Iterate through the child sequences and add their variable root values if
    !the total sequence sums to the target.
    cursor = pstart
    start = pstart
    if (.not. allocated(sequences)) then
       !Sometimes we end up with a sequence that has no kids ever; for example when
       !the first variable is assigned the maximum exponent (equal to powersum). Then
       !we just want to return a single sequence with all other variables set to zero.
       if (this%kidcount .eq. 0) then
          allocate(sequences(1, this%varcount))
       else
          allocate(sequences(this%kidcount, this%varcount))
       end if
       sequences = 0
    end if

    if (associated(this%kids)) then
       !During the initializer, if there are no kids, the kids array could still be initialized
       !to have size zero; this handles that case and sets the first variable in the sequence to
       !have the root. See the comment above for static allocation of sequences(1, ...)
       if (size(this%kids, 1) .eq. 0) then
          sequences(cursor, varindex) = this%root
       else
           do j = 1, size(this%kids,1)
              !Here is where the recursion happens; we add the sequence of this variable's
              !children to the right of this root.
              call this%kids(j)%expand(sequences, cursor, varindex+1)
              cursor = cursor + this%kids(j)%kidcount
              if (varindex .eq. this%varcount-1) cursor = cursor+1
              do k = start, cursor-1
                 sequences(k,varindex) = this%root
              end do
              start = cursor
           end do
        end if
    else
       sequences(cursor,varindex) = this%root
    end if
  end subroutine expand

  !!<summary>Initializes the empty product of multinomials.</summary>
  !!<parameter name="coefficient" regular="true">The coefficient multiplying the product 
  !!of multinomial.</parameter>
  !!<parameter name="targets" regular="true">An array of specific exponents on the 
  !!variables that will contribute to the number of unique colorings.</parameter>
  !!<parameter name="multinomials" regular="true">A mx2 array defining the power and exponent 
  !!of each multinomial in the product. Each row is a multinomial definition. DEALLOCATES 
  !!this variable.</parameter>
  subroutine product_initialize(this, coefficient, targets, multinomials)
    class(Product) :: this
    integer, intent(in) :: coefficient
    integer, allocatable, intent(in) :: targets(:)
    integer, allocatable, intent(inout) :: multinomials(:,:)

    integer :: i

    this%coefficient = coefficient
    this%targets = targets
    allocate(this%multinoms(size(multinomials,1)))
    do i=1, size(multinomials,1)
       call this%multinoms(i)%initialize(multinomials(i,1), multinomials(i,2))
    end do

    !Change the multinomial array to be embedded inside of the product class.
    call move_alloc(multinomials, this%multinom_rd)
  end subroutine product_initialize

  !!<summary>Returns the coefficient of the term with the target exponents if all the multinomials
  !!in the product were expanded and had their terms collected.</summary>
  integer(lint) function coeff(this)
    class(Product), intent(in) :: this

    !!<local name="possibles">
    !!  <dimension index="1">Multinomials in the product.</dimension>
    !!  <dimension index="2">Possible exponents for each of the variables in the multinomial.</dimension>
    !!  <summary>2D array of possible exponents.</summary>
    !!</local>
    !!<local name="mnseq">An array of possible sequences for each multinomial. The 'mnseq' index
    !!runs over the multinomials in the product. The 2D ragged-array in the instance has sequences
    !!as rows and possible variable exponents as columns.</local>
    !!<parameter name="varseq">The algorithm first looks at possible combinations of exponents on
    !!the *first* variable. We then need to construct possible sequences for the remaining variables
    !!in each multinomial. This variable is a temporary Sequence() instance for a single multinomial
    !!which generates an expanded list of possible sequences for that multinomial. It gets cleaned
    !!up after each multinomial.</parameter>
    class(vararray), allocatable :: possibles(:)
    class(vararray2d), allocatable, target :: mnseq(:)
    class(Sequence), pointer :: varseq

    integer :: i, j
    !!<local name="seq0">Array of possible exponents on the *first* variable only. Each row
    !!is a possible set of exponents; each column has the exponent on the *first* variable
    !!for each multinomial in the product; i.e. there are size(this%multinoms) columns.</local>
    !!<local name="expandedseq">Array of possible variable exponents for a *single* multinomial
    !!in the product. Rows are possible sequences of exponents; columns have the exponent for each
    !!variable in the polynomial; i.e. there are size(this%targets) columns.</local>
    !!<local name="seq">Represents a single row in the 'seq0' array, which gets dealt with
    !!one at a time.</local>
    integer, allocatable :: seq0(:,:), expandedseq(:,:)
    integer :: seq(size(this%multinoms, 1))

    allocate(possibles(size(this%multinoms,1)))
    allocate(mnseq(size(this%multinoms, 1)))

    !If this is an isolated multinomial, we only need to check the coefficient of the target term.
    if (size(this%multinoms,1) == 1) then
       if (all(this%multinoms(1)%power-this%targets .gt. 0)) then
          coeff = 0          
       else 
          coeff = this%multinoms(1)%nchoosekm(this%targets)*this%coefficient 
       end if           
    end if

    !Get a list of the possible exponents for each variable in each of the multinomials.
    !We start with the first variable and choose only those combinations of exponents
    !across *all* the multinomials that give the correct target exponent for that variable.
    do i = 1, size(this%multinoms,1)
       call possibles(i)%init(this%multinoms(i)%possible_powers)
    end do

    call cproduct(possibles, seq0)
    !Next, we construct Sequence instances for each of the first variable compatible
    !possibilities and follow them through to the other variables.
    coeff = 0
    allocate(varseq)

    do i = 1, size(seq0,1)
       !Each sequence calculated from the first variable has an entry for each multinomial
       !in this product. The Sequence instances construct smart sequences for the remaining
       !variables in each multinomial separately.
       seq = seq0(i,:)
       do j = 1, size(seq,1)
          call varseq%initialize(seq(j), possibles(j)%items, 2, this%multinoms(j)%powersum, this%targets)
          call varseq%expand(expandedseq, 1, 1)
          call mnseq(j)%init(expandedseq, alloc=.true.)
          !Clean up the varsequence and expanded sequence variable.
          call varseq%finalize()
          deallocate(expandedseq)
       end do

       coeff = coeff + this%sum_sequence(mnseq)
       do j = 1, size(seq,1)    
          call mnseq(j)%finalize()
       end do
    end do

    coeff = coeff*this%coefficient
  end function coeff

  !!<summary>Sums all the possible combinations of relevant sequences based of the variable sequence
  !!lists specified.</summary>
  !!<parameter name="mnseq">a list of possible variable sequences in each 
  !!multinomial (one for each multinomial) that might contribute to the correct target variable.</parameter>
  integer(lint) function sum_sequence(this, mnseq)
    class(Product), intent(in) :: this
    class(vararray2d), intent(in) :: mnseq(size(this%multinoms,1))

    !!<local name="passin">Array of sequence sets (one for each multinomial) that will have their
    !!cartesian product taken.</local>
    class(vararray), allocatable :: passin(:)
    
    integer :: i, j, z
    !!<local name="seqprod">Temporary variable for computing the product over adjacent term coefficients.</local>
    !!<local name="allseqs">Cartesian product of the possible sequences from each multinomial in the product.</local>
    !!<local name="seqlj">Re-constructed S_lj for a single result in 'allseqs' using the sequences 
    !!from 'mnseq'.</local>
    !!<local name="expsum">Temporary sum of the exponents across all adjacent sequences in the product 
    !!(per variable).</local>
    !!<local name="irange">Temporary range from 1..len(mnesq) for constructing the cartesian product.</local>
    integer(lint) :: seqprod
    integer, allocatable :: allseqs(:,:)
    integer :: seqlj(size(this%multinoms, 1), size(this%targets, 1)), expsum(size(this%targets, 1))
    integer, allocatable :: irange(:)

    sum_sequence = 0
    allocate(passin(size(mnseq,1)))
    !We could speed this up by keeping the ranges between loop iterations; if a two of the mnsequences have
    !the same number of elements, we could re-use the range. For now its not worth the extra complexity.
    do i=1, size(this%multinoms,1)
       allocate(irange(mnseq(i)%length))
       irange = (/( j, j=1, mnseq(i)%length, 1 )/)
       call passin(i)%init(irange,alloc=.true.)
       deallocate(irange)
    end do

    call cproduct(passin, allseqs)
    !Each row in allseqs now has a list of integer *indices* from mnseq that need to be used for the
    !coefficient calculations.
    do i = 1, size(allseqs, 1)
       !First we loop over the multinomials and extract the sequence corresponding to the index determined
       !by allseqs.
       do j = 1, size(this%multinoms, 1)
          seqlj(j,:) = mnseq(j)%items(allseqs(i,j),:)
       end do
       !Now that we have the sequences, we sum them up one variable at a time to make sure that they produce the
       !correct target term.
       do j=1, size(this%targets, 1)
          expsum(j) = sum(seqlj(:,j))
       end do
       if (all(expsum .eq. this%targets)) then
          seqprod = 1
          do z = 1, size(seqlj, 1)
             seqprod = seqprod*this%multinoms(z)%nchoosekm(seqlj(z,:))
          end do
          sum_sequence = sum_sequence + seqprod
       end if
    end do
  end function sum_sequence

  !!<summary>Sets up the multinomial.</summary>
  !!<parameter name="power">The value of r for the r-cycle. Also the power on each of the variables
  !!in the multinomial.</parameter>
  !!<parameter name="exponent">The exponent that the entire multinomial is raised to.</parameter>
  subroutine multinomial_initialize(this, power, exponent)
    class(Multinomial) :: this
    integer, intent(in) :: exponent, power
    integer :: j

    this%power = power
    this%exponent = exponent
    this%powersum = power*exponent
    this%possible_powers = (/(j,j=0,power*exponent,power)/)
  end subroutine multinomial_initialize

  !!<summary>Returns the number of different ways to partition an n-element
  !!set into disjoint subsets of sizes k1, ..., km.</summary>
  !!<parameter name="sequencei">an un-normed tuple of form (k1, k2, k3).</parameter>
  integer(lint) function nchoosekm(this, sequencei)
    class(Multinomial), intent(in) :: this
    integer, intent(in) :: sequencei(:)
    integer :: j
    integer :: nsum
    integer :: normseq(size(sequencei,1))

    if (.not. all(mod(sequencei, this%power) == 0)) then
       nchoosekm = 0
    else
       normseq = sequencei/this%power
       nchoosekm = 1
       do j = 1, size(sequencei,1)
          nsum = sum(normseq(1:j))
          nchoosekm = nchoosekm*nchoosek(nsum, normseq(j))
       end do
    end if
  end function nchoosekm

  !!<summary>This implementation was taken from "Binomial Coefﬁcient Computation: Recursion 
  !!or Iteration?" by Yannis Manolopoulos, ACM SIGCSE Bulletin InRoads, Vol.34, No.4, 
  !!December 2002. http://delab.csd.auth.gr/papers/SBI02m.pdf It is supposed to be robust 
  !!against large, intermediate values and to have optimal complexity.</summary>
  !!<parameter name="n" regular="true"></parameter>
  !!<parameter name="k" regular="true"></parameter>
  function nchoosek(n,k)
    integer(lint) :: nchoosek
    integer, intent(in) :: n, k

    !!<local name="t">Variable to accumulate answer in.</local>
    integer(lint) :: t
    integer :: i

    if (k < 0 .or. k > n) then
       nchoosek = 0
       return
    end if
    if (k==0 .and. n == 0) then
       nchoosek = 1
       return
    end if
    t = 1
    if (k<n-k) then
       do i = n, n-k+1, -1
          t = t*i/(n-i+1)
       end do
    else
       do i = n, k+1, -1
          t = t*i/(n-i+1)
       end do
    end if
    nchoosek = t
  end function nchoosek

  !!<summary>Uses a group and concentrations to find the number of unique arrangements as 
  !!described by polya.</summary>
  !!<parameter name="concentrations" regular="true">A list of integers specifying how many of 
  !!each coloring should be present in each of the enumerated lists.</parameter>
  !!<parameter name="group" regular="true">Group operations for permuting the colorings. Rows are operations,
  !!columns are the integer permutations describing the operation.</parameter>
  !!<parameter name="polynomials">A list of integers where the index is the power on the 
  !!variables in the multinomial and the value is the degeneracy of that multinomial in 
  !!the group operation.</parameter>
  !!<parameter name="decompose" regular="true">Specifies whether the r-cycle structure of the
  !!polynomials array should be found explicitly.</parameter>
  integer(lint) function polya(concentrations, group, polynomials, decompose)
    integer, allocatable, intent(in) :: concentrations(:)
    integer, pointer, intent(in) :: group(:,:)
    integer, optional, intent(inout) :: polynomials(size(group,1), size(group,2))
    logical, optional :: decompose

    !!<local name="polyndict">A dictionary of the unique products of multinomials in the polya calculation.</local>
    !!<local name="mult">A mx2 array of the powers and degeneracies of each multinomial in a product.</local>
    !!<local name="operation">Temporary variable for holding a single group operation.</local>
    !!<local name="visited">Temporary array for keeping track of which operation elements have already been
    !!visited when trying to determine the r-cycle structure of the group operation.</local>
    class(Product), allocatable :: polyndict(:)
    integer, allocatable :: mult(:,:)
    integer :: operation(size(group,2)), visited(size(group,2))
    !!<local name="cursor">Points to the current element in the group operation that is being visited.</local>
    !!<local name="vindex">Temporary, used to find the location of the next operation to process when finding
    !!the r-cycles.</local>
    !!<local name="powers">Temporary, the value r for the r-cycle.</local>
    !!<local name="m">The number of unique multinomials in the product.</local>
    !!<local name="mindex">The index of the degeneracy for the lowest valued r-cycle left in the polynomials
    !!array. Corresponds to the value r.</local>
    !!<local name="dupindex">The index of the Product instance in polyndict that has the same multinomial set
    !!as the one for the current group operation.</local>
    !!<parameter name="pused">The number of unique Product instances created so far.</parameter>
    integer :: cursor(1), vindex(1), powers
    integer :: i, j, m, mindex, dupindex, pused = 0
    logical :: ldecomp

    !This checks that the specified concentration matches the number of sites the group is 
    !acting on.
    if (sum(concentrations) .ne. size(group,2)) then
       stop "Error, the concentrations don't add up to the number of sites being used!"
    end if

    !These checks allow a pre-compiled r-cycle structure for the group to be passed in
    !as the form of the polynomials array.
    if (.not. present(polynomials)) then
       polynomials = 0
       ldecomp = .true.
    else
       if (present(decompose)) then
          ldecomp = decompose
          polynomials = 0
       else
          ldecomp = .false.
       end if
    end if

    allocate(polyndict(size(group,1)))
    do i=1, size(group,1)
       if (.not. present(polynomials) .or. ldecomp) then
          operation = group(i,:)
          visited = 0
          do while(any(visited == 0))
             cursor = minloc(visited)
             vindex = minloc(visited)
             powers = 1
             visited(cursor(1)) = 1
             cursor = operation(cursor(1))
             do while(cursor(1) /= vindex(1))
                visited(cursor(1)) = 1
                powers = powers + 1
                cursor = operation(cursor(1))
             end do
             
             polynomials(i,powers) = polynomials(i,powers) + 1          
          end do
       end if

       !We now know the powers and degeneracies of all the r-cycles in the group operation.
       !Create the 2D multinomial array that defines r,d for each M_j^r.
       m = count(polynomials(i,:) .gt. 0)
       allocate(mult(m,2))
       mindex = 1
       do j=1, size(polynomials, 2)
          if (polynomials(i,j) .gt. 0) then
             mult(mindex,1) = j
             mult(mindex,2) = polynomials(i,j)
             mindex = mindex + 1
          end if
       end do

       !Find out if any of the existing products has this multinomial set already.
       dupindex = 0
       print *, "pused ", pused
       do j=1, pused
          if (all(shape(polyndict(j)%multinom_rd) .eq. shape(mult))) then
             print *, "j ", j
             print *, "s1", shape(polyndict(j)%multinom_rd)
             print *, "s2", shape(mult)
             print *, "equal?", (shape(polyndict(j)%multinom_rd) .eq. shape(mult))
             print *, "c1", polyndict(j)%multinom_rd
             print *, "c2", mult
             print *, "equalc?",(polyndict(j)%multinom_rd .eq. mult)
             if (all(polyndict(j)%multinom_rd .eq. mult)) then
                polyndict(j)%coefficient = polyndict(j)%coefficient + 1
                dupindex = j
             end if
          end if
       end do

       if (dupindex == 0) then
          !We found no duplicate product in the set already, just initialize the next empty slot
          pused = pused + 1
          call polyndict(pused)%initialize(1, concentrations, mult)
       end if
       if (allocated(mult)) then 
          deallocate(mult)
       end if
    end do

    polya = 0
    do i = 1, pused
      polya = polyndict(i)%coeff() + polya
    end do    
    polya = polya/size(group,1)
  end function polya
end MODULE classes
