Module HNF_profiler
use num_types
use derivative_structure_generator, only: get_all_HNFs, remove_duplicate_lattices, &
     & get_dvector_permutations
use vector_matrix_utilities, only: minkowski_reduce_basis, norm, matrix_inverse
use enumeration_types
implicit none
public get_HNFs

CONTAINS

  !!<summary>Lists the HNFs for the given symmetry group.</summary>
  !!<parameter name="A" regular="true">The parent lattice for the
  !!system.</parameter>
  !!<parameter name="n_" regular="true">The determinant of the desired
  !!HNF.</parameter>
  !!<parameter name="eps_" regular="true">The floating point
  !!tolerance.</parameter>
  !!<parameter name="hnfs" regular="true">The returned list of
  !!HNFs.</parameter>
  !!<parameter name="r_mins" regular="true">The r_min values for the
  !!HNFs.</parameter>
  subroutine get_HNFs(A,hnfs,r_mins,n_,eps_)
    real(dp), intent(in) :: A(3,3)
    integer, intent(out) :: hnfs(100,3,3)
    real(dp), intent(out) :: r_mins(100)
    integer, intent(in), optional :: n_
    real(dp), intent(in), optional :: eps_

    integer :: n_min, n_max, LatDim, n_hnfs, i, j, count_i, count_j
    real(dp) :: eps, max_rmin, this_rmin
    type(RotPermList) :: dRPList
    type(RotPermList), pointer :: RPlist(:)
    integer, pointer :: all_hnfs(:,:,:), unq_hnfs(:,:,:), degeneracy_list(:)
    type(opList), pointer :: fixop(:)
    real(dp), pointer :: latts(:,:,:), d(:,:)
    real(dp) :: reduced_latt(3,3), norms(3), inv_lat(3,3)
    integer :: temp_hnfs(100,3,3)

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

    call get_dvector_permutations(A,d,dRPList,LatDim,eps)
    call matrix_inverse(A,inv_lat)

    do i=n_min,n_max
       call get_all_HNFs(i,all_hnfs)
       call remove_duplicate_lattices(all_hnfs,LatDim,A,d,dRPList,unq_hnfs,fixop, RPList, latts, degeneracy_list,eps)
       count_i = 0
       max_rmin = 0
       temp_hnfs = 0
       do j=1,size(latts,3)
          call minkowski_reduce_basis(latts(:,:,j),reduced_latt,eps)
          norms = norm(reduced_latt)
          this_rmin = minval(norms)
          if (count_i ==0) then
             count_i = count_i + 1
             max_rmin = this_rmin
             temp_hnfs(count_i,:,:) = matmul(inv_lat,reduced_latt)
          else if (abs(this_rmin-max_rmin)<eps) then
             count_i = count_i + 1
             temp_hnfs(count_i,:,:) = NINT(matmul(inv_lat,reduced_latt))
          else if (this_rmin > max_rmin) then
             count_i = 1
             temp_hnfs = 0
             max_rmin = this_rmin
             temp_hnfs(count_i,:,:) = matmul(inv_lat,reduced_latt)
          end if
       end do

       do j=1,count_i
          n_hnfs = n_hnfs + 1
          hnfs(n_hnfs,:,:) = temp_hnfs(j,:,:)
          r_mins(n_hnfs) = max_rmin
          if (n_hnfs==100) exit 
       end do
       deallocate(latts,unq_hnfs,RPList,fixOp)
       if (n_hnfs==100) exit 
    end do
  end subroutine get_HNFs

end Module HNF_profiler
