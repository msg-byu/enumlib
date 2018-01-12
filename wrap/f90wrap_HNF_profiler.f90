! Module hnf_profiler defined in file ../aux_src/HNF_profiler.f90

subroutine f90wrap_get_hnfs(a, hnfs, r_mins, n_, eps_)
    use hnf_profiler, only: get_hnfs
    implicit none
    
    real(8), intent(in), dimension(3,3) :: a
    integer, dimension(100,3,3), intent(inout) :: hnfs
    real(8), dimension(100), intent(inout) :: r_mins
    integer, intent(in), optional :: n_
    real(8), intent(in), optional :: eps_

    call get_hnfs(A=a, hnfs=hnfs, r_mins=r_mins, n_=n_, eps_=eps_)
end subroutine f90wrap_get_hnfs

! End of module hnf_profiler defined in file ../aux_src/HNF_profiler.f90

