! Module hnf_profiler defined in file ../aux_src/HNF_profiler.f90

subroutine f90wrap_get_hnfs(a, basis, hnfs, r_mins, n_, eps_, n0, n1)
    use hnf_profiler, only: get_hnfs
    implicit none
    
    real(8), intent(in), dimension(3,3) :: a
    real(8), intent(in), dimension(n0,n1) :: basis
    integer, dimension(3,3,100), intent(inout) :: hnfs
    real(8), dimension(100), intent(inout) :: r_mins
    integer, intent(in), optional :: n_
    real(8), intent(in), optional :: eps_
    integer :: n0
    !f2py intent(hide), depend(basis) :: n0 = shape(basis,0)
    integer :: n1
    !f2py intent(hide), depend(basis) :: n1 = shape(basis,1)
    call get_hnfs(A=a, basis=basis, hnfs=hnfs, r_mins=r_mins, n_=n_, eps_=eps_)
end subroutine f90wrap_get_hnfs

! End of module hnf_profiler defined in file ../aux_src/HNF_profiler.f90

