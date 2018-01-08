import _supercell
import f90wrap.runtime
import logging

class Hnf_Profiler(f90wrap.runtime.FortranModule):
    """
    Module hnf_profiler
    
    
    Defined at ../aux_src/HNF_profiler.f90 lines 1-97
    
    """
    @staticmethod
    def get_hnfs(a, hnfs, r_mins, n_=None, eps_=None):
        """
        get_hnfs(a, basis, hnfs, r_mins[, n_, eps_])
        
        
        Defined at ../aux_src/HNF_profiler.f90 lines 24-97
        
        Parameters
        ----------
        a : float array
        basis : float array
        hnfs : int array
        r_mins : float array
        n_ : int
        eps_ : float
        
        """
        _supercell.f90wrap_get_hnfs(a=a, hnfs=hnfs, r_mins=r_mins, n_=n_, \
            eps_=eps_)
    
    _dt_array_initialisers = []
    

hnf_profiler = Hnf_Profiler()

