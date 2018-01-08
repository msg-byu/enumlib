import _supercell
import f90wrap.runtime
import logging
import numpy as np

def get_supers(atoms,n_=None,eps_=None):
    """Generates the hnfs for the given atoms object from n=2 multiples to
    n_ (default 100) multiples of the parent cell.

    Args:
        atoms (quipy.atoms.Atoms): an atoms object that the supercell 
            matrices will be found for.
        n_ (int): optional maximum number of multiples to find.
        eps_ (float): optional floating point tolerance.
    
    Returns:
        (hnfs,rmin) the integer matrices that create unique supercells of the parent,
            and the rmin_values that match.
    """

    hnfs = np.zeros((100,3,3), dtype=np.int32, order='F')
    rmin = np.zeros(100, dtype=float, order='F')
    Hnf_Profiler.get_hnfs(np.transpose(atoms.cell),hnfs,rmin,n_=n_,eps_=eps_)
    keep = np.where(rmin != 0)[0]

    return (hnfs[keep],rmin[keep])

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

