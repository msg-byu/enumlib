Module HNF_profiler
use num_types
use derivative_structure_generator, only: get_all_HNFs, remove_duplicate_lattices, &
     & get_dvector_permutations
use vector_matrix_utilities, only: minkowski_reduce_basis, norm, matrix_inverse, determinant
use symmetry, only: make_primitive, get_lattice_pointGroup
use enumeration_types
implicit none
public get_HNFs

CONTAINS

  !!<summary>Lists the HNFs for the given symmetry group.</summary>
  !!<parameter name="A" regular="true">The parent lattice for the
  !!system.</parameter>
  !!<parameter name="atom_pos" regular="true">The atomic basis in
  !!cartesian coordinates.</parameter>
  !!<parameter name="atom_types" regular="true">The atomic occupansies
  !!for each atom in the cell (atom A = 0, B = 1 ....).</parameter>
  !!<parameter name="n_" regular="true">The determinant of the desired
  !!HNF.</parameter>
  !!<parameter name="eps_" regular="true">The floating point
  !!tolerance.</parameter>
  !!<parameter name="hnfs" regular="true">The returned list of
  !!HNFs.</parameter>
  !!<parameter name="r_mins" regular="true">The r_min values for the
  !!HNFs.</parameter>
  !!<parameter name="r_maxs" regular="true">The r_max values for the
  !!HNFs.</parameter>
  !!<parameter name="dets" regular="true">The determinants of the
  !!HNFs returned.</parameter>
  !!<parameter name="pgs" regular="true">The size of the
  !!point-group for the superlattice.</parameter>
  subroutine get_HNFs(A, atom_pos, atom_types, hnfs, r_mins, r_maxs, pgs, dets, n_, eps_)
    real(dp), intent(in) :: A(3,3), atom_pos(:,:)
    integer, intent(in) :: atom_types(:)
    integer, intent(out) :: hnfs(100,3,3)
    real(dp), intent(out) :: r_mins(100), r_maxs(100)
    integer, intent(out) :: pgs(100), dets(100)
    integer, intent(in), optional :: n_
    real(dp), intent(in), optional :: eps_

    integer :: n_min, n_max, LatDim, n_hnfs, i, j, count_i
    real(dp) :: eps, max_rmin, this_rmin, this_rmax
    type(RotPermList) :: dRPList
    type(RotPermList), pointer :: RPlist(:)
    integer, pointer :: all_hnfs(:,:,:), unq_hnfs(:,:,:), degeneracy_list(:)
    type(opList), pointer :: fixop(:)
    real(dp), pointer :: latts(:,:,:), d(:,:), pg_Ops(:,:,:), atomPos(:,:)
    real(dp) :: reduced_latt(3,3), norms(3), inv_lat(3,3), prim_lat(3,3)
    integer :: temp_hnfs(100,3,3), det_fix, js(100), HNF(3,3),temp_mat(3,3)
    integer, pointer :: atomTypes(:)
    real(dp), allocatable, target :: temp_pos(:,:)
    integer, allocatable, target :: temp_types(:)

    ! It only makes sense to use this code for bulk crystals and not
    ! surfaces so we hardcode the LatDim to be 3.
    LatDim = 3
    n_min = 2
    if (.not. present(n_)) then
       n_max = 100
    else
       n_max = n_
    end if

    if (.not. present(eps_)) then
       eps = 1E-3
    else
       eps = eps_
    end if

    allocate(d(3,1))
    d = 0.0_dp
    hnfs = 0
    n_hnfs = 0
    r_mins = 0

    ! This is so we can fix the determinants of the HNFs for left
    ! handed systems to be positive.
    if (determinant(A) < 0) then
       det_fix = -1
    else
       det_fix = 1
    end if

    prim_lat = A

    ! Here we allocate temp variables for the atom types and positions
    ! and then associate them with pointers that will be passed into
    ! make_primitive. This is necessary because F90wrap cannot have
    ! pointers declared in the subroutine call derictley and so we
    ! have to define the pointers internally. The pointers also can't
    ! be set equal to the input variables becaus it causes segfaults
    ! when the code is run, hence the temporary variables.
    allocate(temp_pos(size(atom_pos,1),size(atom_pos,2)),temp_types(size(atom_types,1)))
    temp_pos = atom_pos
    temp_types = atom_types 
    atomPos => temp_pos
    atomTypes => temp_types
    call make_primitive(prim_lat, atomTypes, atomPos, .False., eps_=eps)

    call get_dvector_permutations(prim_lat, d, dRPList, LatDim, eps)
    call matrix_inverse(prim_lat, inv_lat)

    do i=n_min,n_max
       ! Find all the HFNs then remove those that are equivalent.
       call get_all_HNFs(i, all_hnfs)
       call remove_duplicate_lattices(all_hnfs, LatDim, prim_lat, d, dRPList, unq_hnfs, &
            fixop, RPList, latts, degeneracy_list,eps)
       count_i = 0
       max_rmin = 0
       temp_hnfs = 0
       ! We need to search the HNFs for the one that has the maximum
       ! r_min. For now we keep all ties but we'll need to change this
       ! latter.
       do j=1,size(latts,3)
          call minkowski_reduce_basis(latts(:,:,j), reduced_latt, eps)
          norms = norm(reduced_latt)
          this_rmin = minval(norms)
          this_rmax = maxval(norms)
          if (count_i ==0) then
             count_i = count_i + 1
             max_rmin = this_rmin
             temp_hnfs(count_i,:,:) = NINT(matmul(inv_lat, reduced_latt))
             js(count_i) = j
          else if (abs(this_rmin-max_rmin) < eps) then
             count_i = count_i + 1
             temp_hnfs(count_i,:,:) = NINT(matmul(inv_lat, reduced_latt))
             js(count_i) = j
          else if (this_rmin > max_rmin) then
             js = 0
             count_i = 1
             temp_hnfs = 0
             max_rmin = this_rmin
             temp_hnfs(count_i,:,:) = NINT(matmul(inv_lat, reduced_latt))
             js(count_i) = j
          end if
       end do
       
       do j=1,count_i
          ! We repeat the minkowski reduction here to get the rmax
          ! value and retrieve the point group for the lattice. Later
          ! this will be moved into the hnf selection above. We may
          ! also switch to finding the spaceGroup instead of the point
          ! group at later down the road, for now the lattice point
          ! group gives us a good heuristic.
          call minkowski_reduce_basis(latts(:,:,js(j)), reduced_latt, eps)
          call get_lattice_pointGroup(reduced_latt, pg_Ops)
          norms = norm(reduced_latt)
          this_rmax = maxval(norms)
          n_hnfs = n_hnfs + 1
          hnfs(n_hnfs,:,:) = det_fix*temp_hnfs(j,:,:)
          r_mins(n_hnfs) = max_rmin
          r_maxs(n_hnfs) = this_rmax
          dets(n_hnfs) = i
          pgs(n_hnfs) = size(pg_Ops,3)
          if (n_hnfs==100) exit 
       end do
       deallocate(latts, unq_hnfs, RPList, fixOp)
       js = 0
       if (n_hnfs==100) exit 
    end do
  end subroutine get_HNFs

end Module HNF_profiler
