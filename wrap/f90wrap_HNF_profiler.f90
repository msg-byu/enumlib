! Module hnf_profiler defined in file ../aux_src/HNF_profiler.f90

subroutine f90wrap_get_hnfs(a, atom_pos, atom_types, hnfs, r_mins, r_maxs, pgs, &
    dets, n0, n1, n2, n_, eps_)
    use hnf_profiler, only: get_hnfs
    implicit none

    real(8), intent(in), dimension(3,3) :: a
    real(8), intent(in), dimension(n0,n1), allocatable :: atom_pos(:,:)
    integer, intent(inout), dimension(n2) :: atom_types
    integer, dimension(100,3,3), intent(inout) :: hnfs
    real(8), dimension(100), intent(inout) :: r_mins
    real(8), dimension(100), intent(inout) :: r_maxs
    integer, dimension(100), intent(inout) :: pgs
    integer, dimension(100), intent(inout) :: dets
    integer, intent(in), optional :: n_
    real(8), intent(in), optional :: eps_
    integer :: n0
    !f2py intent(hide), depend(atom_pos) :: n0 = shape(atom_pos,0)
    integer :: n1
    !f2py intent(hide), depend(atom_pos) :: n1 = shape(atom_pos,1)
    integer :: n2
    !f2py intent(hide), depend(atom_types) :: n2 = shape(atom_types,0)
    call get_hnfs(A=a, atom_pos=atom_pos, atom_types=atom_types, hnfs=hnfs, &
        r_mins=r_mins, r_maxs=r_maxs, pgs=pgs, dets=dets, n_=n_, eps_=eps_)
end subroutine f90wrap_get_hnfs

! End of module hnf_profiler defined in file ../aux_src/HNF_profiler.f90
