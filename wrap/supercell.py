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

    # For now we want to keep track of all these outputs, later we may
    # trim down the list.
    hnfs = np.zeros((100,3,3), dtype=np.int32, order='F')
    rmin = np.zeros(100, dtype=float, order='F')
    rmax = np.zeros(100, dtype=float, order='F')
    pgs = np.zeros(100, dtype=np.int32, order='F')
    dets = np.zeros(100, dtype=np.int32, order='F')

    atoms_nums = atoms.numbers
    u_atoms = np.unique(atoms_nums)
    atom_types = []
    for i in atoms_nums:
        atom_types.append(np.where(u_atoms==i)[0][0])    
    
    Hnf_Profiler.get_hnfs(np.transpose(atoms.cell), np.transpose(atoms.positions),
                          atom_types, hnfs, rmin, rmax, pgs, dets, n_=n_, eps_=eps_)
    # We have more entries than needed. We only want to keep the non-zero ones.
    keep = np.where(rmin != 0)[0]

    # We apply a round to trim down the outputs for rmin and rmax.
    return (hnfs[keep],np.round(rmin[keep],3),np.round(rmax[keep],3),pgs[keep],dets[keep])

class Hnf_Profiler(f90wrap.runtime.FortranModule):
    """
    Module hnf_profiler
    
    
    Defined at ../aux_src/HNF_profiler.f90 lines 1-141
    
    """
    @staticmethod
    def get_hnfs(a, atom_pos, atom_types, hnfs, r_mins, r_maxs, pgs, dets, n_=None, \
        eps_=None):
        """
        get_hnfs(a, atom_pos, atom_types, hnfs, r_mins, r_maxs, pgs, dets[, n_, eps_])
        
        
        Defined at ../aux_src/HNF_profiler.f90 lines 34-141
        
        Parameters
        ----------
        a : float array
        atom_pos : float array
        atom_types : int array
        hnfs : int array
        r_mins : float array
        r_maxs : float array
        pgs : int array
        dets : int array
        n_ : int
        eps_ : float
        
        """
        _supercell.f90wrap_get_hnfs(a=a, atom_pos=atom_pos, atom_types=atom_types, \
            hnfs=hnfs, r_mins=r_mins, r_maxs=r_maxs, pgs=pgs, dets=dets, n_=n_, \
            eps_=eps_)
    
    _dt_array_initialisers = []
    

hnf_profiler = Hnf_Profiler()
