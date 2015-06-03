#if 0
!----------------------------------------------------------
! This is a test driver module meant for testing the derivative_structure_generator.
! 
! This will test the individual functions and give output where code is
! broken for debugging purposes.
!
! Author: Derek Carr
! Email: zaphinath@gmail.com
!
!----------------------------------------------------------
#endif

MODULE test_driver_dsg
use derivative_structure_generator

implicit none


CONTAINS
SUBROUTINE test_gen_multilattice_derivatives()
	!call gen_multilattice_derivatives(title, parLV, nDFull, dFull, k, nMin, nMax, pLatTyp, eps, full,&
    !& labelFull,digitFull,equivalencies,conc_check,conc_ElementN, conc_Range)
END SUBROUTINE test_gen_multilattice_derivatives


SUBROUTINE test_all_dsg()
	write(*,'(A)')	"Case with test call"
END SUBROUTINE test_all_dsg

END MODULE test_driver_dsg


