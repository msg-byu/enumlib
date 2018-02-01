"""This scripts produces the VASP style POSCAR."""
from __future__ import print_function
import warnings
import numpy as np

try:
    from termcolor import cprint
except ImportError:
    cprint = print

# The dictionary of all the elements on the periodic table
element_volume ={"H":37.2958,"He":32.1789,"Li":21.2543,"Be":8.49323,"B":7.24205,"C":5.68741,
                 "N":46.6002,"O":22.2802,"F":17.0258,"Ne":21.7346,"Na":23.2596,"Mg":23.3928,
                 "Al":16.6075,"Si":7.8511,"P":9.1459,"S":17.1672,"Cl":35.2074,"Ar":36.3829,
                 "K":71.5278,"Ca":43.4353,"Sc":25.6478,"Ti":18.1565,"V":13.7718,"Cr":11.9439,
                 "Mn":19.3207,"Fe":11.82,"Co":11.1838,"Ni":10.9036,"Cu":11.7615,"Zn":13.311,
                 "Ga":18.4496,"Ge":19.3638,"As":20.4270,"Se":58.6173,"Br":33.3170,"Kr":46.7873,
                 "Rb":87.3384,"Sr":56.1889,"Y":33.0792,"Zr":23.8327,"Nb":17.9685,"Mo":15.6279,
                 "Tc":14.5458,"Ru":13.9206,"Rh":13.718,"Pd":14.716,"Ag":17.1045,"Cd":18.7161,
                 "In":26.6861,"Sn":29.3238,"Sb":27.1733,"Te":62.3227,"I":24.3807,"Xe":59.582,
                 "Cs":110.723,"Ba":63.253,"Hf":23.1748,"Ta":18.1323,"W":15.7772,"Re":14.8694,
                 "Os":14.5485,"Ir":14.1558,"Pt":15.0591,"Au":16.9793,"Hg":27.6914,"Tl":29.2949,
                 "Pd":30.3218,"Bi":31.2849}

verbosity = None
"""The verbosity level of messages being printed by the module."""
quiet = None
"""When in quiet mode, no text is ever printed to terminal unless its
verbosity level is explicitly specified.
"""
cenum = {
    "cerrs": 0,
    "cwarn": 1,
    "cinfo": 2,
    "cgens": 3,
    "cstds": 4,
    "cokay": 5
}
"""Dict of the various colors available for coloring specific
lines in the arb().
"""
icols = ["red", "yellow", "cyan", "blue", "white", "green"]
nocolor = False
"""When true, the colored outputs all use the regular print() instead 
so that the stdout looks ordinary.
"""


def _common_parser():
    """Returns a parser with common command-line options for all the scripts
    in the fortpy suite.
    """
    import argparse
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-examples", action="store_true",
                        help="See detailed help and examples for this script.")
    parser.add_argument("-verbose", action="store_true",
                        help="See verbose output as the script runs.")
    parser.add_argument('-action', nargs=1, choices=['save','print'], default='print',
                        help="Specify what to do with the output (print or save)")
    parser.add_argument("-debug", action="store_true",
                        help="Print verbose calculation information for debugging.")

    return parser

bparser = _common_parser()
testmode = False
"""bool: when True, the package is operating in unit test mode, which changes
how plotting is handled.
"""

def deprecated(func):
    '''This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used.'''
    def new_func(*args, **kwargs):
        warnings.warn("Call to deprecated function {}.".format(func.__name__),
                      category=DeprecationWarning)
        return func(*args, **kwargs)
    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    new_func.__dict__.update(func.__dict__)
    return new_func

def exhandler(function, parser):
    """If -examples was specified in 'args', the specified function
    is called and the application exits.
    :arg function: the function that prints the examples.
    :arg parser: the initialized instance of the parser that has the
      additional, script-specific parameters.
    """
    args = vars(bparser.parse_known_args()[0])
    if args["examples"]:
        function()
        return
    if args["verbose"]:
        from liveserial.msg import set_verbosity
        set_verbosity(args["verbose"])

    args.update(vars(parser.parse_known_args()[0]))
    return args

def set_testmode(testing):
    """Sets the package testing mode.
    """
    global testmode
    testmode = testing

def RepresentsInt(s):
    """Determines if a string could be an int.
    :arg s: The string to be tested.
    """    
    try: 
        int(s)
        return True
    except ValueError:
        return False
    
def RepresentsFloat(s):
    """Determines if a string could be an float.
    :arg s: The string to be tested.
    """    
    try: 
        float(s)
        return True
    except ValueError:
        return False

def example(script, explain, contents, requirements, output, outputfmt, details):
    """Prints the example help for the script."""
    blank()
    cprint(script.upper(), "yellow")
    cprint(''.join(["=" for i in range(70)]) + '\n', "yellow")

    cprint("DETAILS", "blue")
    std(explain + '\n')

    cprint(requirements, "red")
    cprint(output, "green")
    blank()

    if details != "":
        std(details)
        blank()

    cprint("OUTPUT FORMAT", "blue")
    std(outputfmt)
    blank()

    cprint("EXAMPLES", "blue")
    for i in range(len(contents)):
        pre, code, post = contents[i]
        std("{}) {}".format(i + 1, pre))
        cprint("    " + code, "cyan")
        if post != "":
            std('\n' + post)
        blank()

def printer(text, color=None, **kwargs):
    """Prints using color or standard print() depending on the value
    of 'nocolor'.
    """
    if nocolor:
        # import sys
        # sys.stdout.write(text + "" if ("end" in kwargs and kwargs["end"] == "") else '\n')
        # sys.stdout.flush()
        print(text, **kwargs)
    else:
        if color is None:
            cprint(text, **kwargs)
        else:
            cprint(text, color, **kwargs)

def arb(text, cols, split):
    """Prints a line of text in arbitrary colors specified by the numeric
    values contained in msg.cenum dictionary.
    """
    stext = text if text[-1] != split else text[0:-1]
    words = stext.split(split)
    for i, word in enumerate(words):
        col = icols[cols[i]]
        printer(word, col, end="")
        if i < len(words)-1:
            printer(split, end="")
        else:
            printer(split)

def set_verbosity(level):
    """Sets the modules message verbosity level for *all* messages printed.
    :arg level: a positive integer (>0); higher levels including more detail.
    """
    global verbosity
    verbosity = level

def set_quiet(is_quiet):
    """Sets whether the messaging system is running in quiet mode. Quiet mode
    only prints messages with explicit verbosity specified. If verbosity==None
    and quiet==True, any message with level >= 0 is printed.
    """
    global quiet
    quiet = is_quiet

def will_print(level=1):
    """Returns True if the current global status of messaging would print a
    message using any of the printing functions in this module.
    """
    if level == 1:
        #We only affect printability using the quiet setting.
        return quiet is None or quiet == False
    else:
        return ((isinstance(verbosity, int) and level <= verbosity) or
                (isinstance(verbosity, bool) and verbosity == True))
    
def warn(msg, level=0, prefix=True):
    """Prints the specified message as a warning; prepends "WARNING" to
    the message, so that can be left off.
    """
    if will_print(level):
        printer(("WARNING: " if prefix else "") + msg, "yellow")

def err(msg, level=-1, prefix=True):
    """Prints the specified message as an error; prepends "ERROR" to
    the message, so that can be left off.
    """
    if will_print(level) or verbosity is None:
        printer(("ERROR: " if prefix else "") + msg, "red")

def info(msg, level=1):
    """Prints the specified message as information."""
    if will_print(level):
        printer(msg, "cyan")

def okay(msg, level=1):
    """Prints the specified message as textual progress update."""
    if will_print(level):
        printer(msg, "green")

def gen(msg, level=1):
    """Prints the message as generic output to terminal."""
    if will_print(level):
        printer(msg, "blue")

def blank(n=1, level=2):
    """Prints a blank line to the terminal."""
    if will_print(level):
        for i in range(n):
            print("")

def std(msg, level=1):
    """Prints using the standard print() function."""
    if will_print(level):
        print(msg)
    
def _gaussian_reduce_two_vectors(U,V,eps):
    """This routine takes two vectors (in three-space) and reduces them to
    form a shortest set (Minkowski reduced). The idea is to subtract
    multiples of U from V so that the new V is as close to the origin
    as any lattice point along the line that passes through U in the
    direction of V. The process is repeated until the new vector isn't
    shorter than the other. It's pretty obvious if you do an example
    by hand. Also see 3.1 of Lecture notes in computer science, ISSN
    0302-974, ANTS - VI: algorithmic number theory, 2004, vol. 3076,
    pp. 338-357 ISBN 3-540-22156-5. Fixes made Apr 2012 GLWH (not sure
    if they made a practical difference though)
    :arg U: a vector
    :arg V: a vector
    :arg eps: finite precision tolerance
    """

    from numpy.linalg import norm
    from numpy import dot

    it = 0
    if norm(U) > (norm(V) - eps):
       # Make sure that the {U,V} are listed in ascending order; ||U||<||V||
       temp = U
       U = V
       V = temp # Keep V as the longest vector

    done = False
    it = 1
    while not done:
        if it > 10: # pragma: no cover
            err("gaussian_reduce_two_vectors failed to converge in 10 iterations")
            exit()
        R = [V[i]-int(round(dot(U,V)/dot(U,U)+1E-10))*U[i] for i in range(3)] #Shorten V as much as possible
        V = U # Swap U and V (so U remains the shortest)
        U = R
        if norm(U) >= (norm(V) - eps):
            done = True
        it += 1

    # Make sure that the {U,V} are listed in ascending order on exit; ||U||<||V||
    temp = U
    U = V
    V = temp
    return U, V

def _minkowski_conditions_check(basis,eps):
    """This function checks the minkowski conditions for a 3D lattice
    basis.
    :arg basis: The atomic basis vectors
    :arg eps: finitie precision tolerance
    """
    from numpy import linalg
    
    b1 = basis[0]
    b2 = basis[1]
    b3 = basis[2]

    minkowski_check = True
    if linalg.norm(b1) > (linalg.norm(b2) + eps):
        minkowski_check = False
        err("Minkowski_condition 1 failed: b1 > b2")

    if linalg.norm(b2) > (linalg.norm(b3) + eps):
        minkowski_check = False
        err("Minkowski_condition 2 failed: b2 > b3")

    if linalg.norm(b2) > (linalg.norm([b1[i]+b2[i] for i in range(len(b1))])+eps):
        minkowski_check = False
        err("Minkowski_condition 3 failed: b2 > b1+b2")

    if linalg.norm(b2) > (linalg.norm([b1[i]-b2[i] for i in range(len(b1))])+eps):
        minkowski_check = False
        err("Minkowski_condition 4 failed: b2 > b1-b2")

    if linalg.norm(b3) > (linalg.norm([b1[i]+b3[i] for i in range(len(b1))])+eps):
        minkowski_check = False
        err("Minkowski_condition 5 failed: b3 > b1+b3")

    if linalg.norm(b3) > (linalg.norm([b3[i]-b1[i] for i in range(len(b1))])+eps):
        minkowski_check = False
        err("Minkowski_condition 6 failed: b3 > b3-b1")

    if linalg.norm(b3) > (linalg.norm([b2[i]+b3[i] for i in range(len(b2))])+eps):
        minkowski_check = False
        err("Minkowski_condition 7 failed: b3 > b2+b3")

    if linalg.norm(b3) > (linalg.norm([b3[i]-b2[i] for i in range(len(b2))])+eps):
        minkowski_check = False
        err("Minkowski_condition 8 failed: b3 > b3-b2")

    if linalg.norm(b3) > (linalg.norm([b1[i]+b2[i]+b3[i] for i in range(len(b2))])+eps):
        minkowski_check = False
        err("Minkowski_condition 9 failed: b3 > b1+b2+b3")

    if linalg.norm(b3) > (linalg.norm([b1[i]-b2[i]+b3[i] for i in range(len(b2))])+eps):
        minkowski_check = False
        err("Minkowski_condition 10 failed: b3 > b1-b2+b3")

    if linalg.norm(b3) > (linalg.norm([b1[i]+b2[i]-b3[i] for i in range(len(b2))])+eps):
        minkowski_check = False
        err("Minkowski_condition 11 failed: b3 > b1+b2-b3")

    if linalg.norm(b3) > (linalg.norm([b1[i]-b2[i]-b3[i] for i in range(len(b2))])+eps):
        minkowski_check = False
        err("Minkowski_condition 12 failed: b3 > b1-b2-b3")

    return minkowski_check

def _reduce_C_in_ABC(A,B,C,eps):
    """This routine takes three vectors, A,B,C, defining a lattice, and
    reduces the last one so that it is as close as possible to the
    origin while remaining in an affine plane, which is parallel to
    the A-B plane but which passes through the end of the C
    vector. See Lecture notes in computer science, ISSN 0302-974, ANTS
    - VI : algorithmic number theory, 2004, vol. 3076, pp. 338-357
    ISBN 3-540-22156-5
    :arg A: a vector
    :arg B: a vector
    :arg C: a vector
    :arg eps: finite precision tolerance
    """
    from numpy import cross, linalg, dot, allclose, matmul, array
    from copy import deepcopy
    from math import floor
    
    oldABC = deepcopy([A,B,C])
    
    # Use Gaussian reduction to reduce the A,B 2D basis so that it is
    # itself Minkowski reduced. If this is done, then the closest
    # lattice point (in A,B plane) to the projection of C (into the
    # A,B plane) is guaranteed to be one of the corners of the unit
    # cell enclosing the projection of C
    (A,B) = _gaussian_reduce_two_vectors(A,B,eps)

    # First thing to do is find the (real, not lattice) point in the
    # affine plane A,B + C that is nearest the origin. Call this T.
    cpdAB = [i/linalg.norm(cross(A,B)) for i in cross(A,B)]
    T = [C[i] - cpdAB[i]*dot(C,cpdAB) for i in range(3)]

    if not allclose(dot(T,cross(A,B)),0,atol=eps,rtol=eps): #pragma: no cover
        err("{} Projection of C into A,B plane failed".format(str(dot(T,cross(A,B)))))

    # Now find the four points of the A,B lattice, in the affine
    # plane, that enclose the point T
    ABC = [A,B,C]
    ABCinv = linalg.inv(ABC)

    LC = [int(floor(i +eps)) for i in matmul(T,ABCinv).tolist()]
    # Compute the distance from T to each of the four corners of the cell and pick
    # the one that is the closest.
    corners = array([[0,0,0],[1,0,0],[0,1,0],[1,1,0]])
    dist = []
    for i in range(0,4):
        temp1 = corners[i] + array(LC)
        temp2 = array(T) -matmul((corners[i] + array(LC)),ABC)
        dist.append(linalg.norm(array(T) -matmul((corners[i] + array(LC)),ABC)))

    idx = dist.index(min(dist))

    if idx == 0:
        temp1 = [corners[0][i] + LC[i] for i in range(3)]
        temp2 = matmul(temp1,ABC).tolist()
        C = [C[i] - temp2[i] for i in range(len(C))]
    elif idx == 1:
        temp1 = [corners[1][i] + LC[i] for i in range(3)]
        temp2 = matmul(temp1,ABC).tolist()
        C = [C[i] - temp2[i] for i in range(len(C))]
    elif idx == 2:
        temp1 = [corners[2][i] + LC[i] for i in range(3)]
        temp2 = matmul(temp1,ABC).tolist()
        C = [C[i] - temp2[i] for i in range(len(C))]
    elif idx == 3:
        temp1 = [corners[3][i] + LC[i] for i in range(3)]
        temp2 = matmul(temp1,ABC).tolist()
        C = [C[i] - temp2[i] for i in range(len(C))]
    else: #pragma: no cover
        err("Case failed in reduce_C_in_ABC"
            "Lattice coordinates in the A,B plane: ".format(' '.join([str(i) for i in LC])))

    ABC = [A,B,C]
    ABCinv = linalg.inv(ABC)
    temp = matmul(list(map(list,zip(*ABCinv))),list(map(list,zip(*oldABC)))).tolist()
    for i in range(3):
        for j in range(3):
            if abs(temp[i][j] - int(round(temp[i][j]))) > eps: #pragma: no cover
                err("Lattice was not preserved in reduce_C_in_ABC")
                exit()

    return A, B, C

def _get_lat_param_element(lat_vecs,n_basis_atoms,element):
    """Finds the lattice parameter of an element for the given set of lattice.
    Args:
        lat_vecs (list): The lattice vectors for the system.
        n_basis_atoms (int): The number of atoms in the atomic basis.
        element (str): The symbol for the element.
    
    Returns:
        lat_param (float): The lattice parameter.
    """

    lat_vol = abs(np.linalg.det(lat_vecs)/n_basis_atoms)
    atom_vol = element_volume[element]/lat_vol

    return atom_vol**(1./3.)   

def _get_lattice_parameter(elements, concentrations, lat_vecs, n_basis_atoms, default_title,remove_zeros=False):
    """Finds the lattice parameters for the provided atomic species using Vagars law.
    Args:
        elements (list of str): A list of the elements in the system.
        default_title (str): The default system title.
        concentrations (list of int): The concentrations of each element.
        lat_vecs (list): The lattice vectors for the system.
        n_basis_atoms (int): The number of atoms in the atomic basis.
    Returns:
        lat_param (float): The lattice parameter of the system.
    Raises:
        ValueError: if the number of elements doesn't match the len of the concentration array.
    """

    if elements is None:
        lat_param = 1.0
        title = default_title
    else:
        if len(elements) != len(concentrations):
            raise ValueError("You have provided {0} element names when {1} elements are present "
                "in the system. Please provide the correct number of elements."
                .format(len(elements),len(concentrations)))

        else:
            title = ""
            lat_param = 0
            for i, elem in enumerate(elements):
                lat_param += concentrations[i]*_get_lat_param_element(lat_vecs,n_basis_atoms,elem)
                if concentrations[i] > 0:
                    title += " {0} ".format(elem)
                elif not remove_zeros:
                    title += " {0} ".format(elem)
            lat_param = float(lat_param) / sum(concentrations)
            title = "{0} {1}\n".format(title,default_title.strip())
    return lat_param, title

def _cartesian2direct(sLV,aBas, eps):
    """This routine takes three lattice vectors and a list of atomic basis
    vector in Cartesian coordinates and converts them to direct
    ("lattice") coordinates.
    :arg sLV: Superlattice vectors (Cartesian coordinates)
    :arg aBas: Atomic positions (cartesian coordinates first, then direct)
    :arg eps: Finite precision tolerance.
    """
    from numpy import linalg, matmul, array
    
    nAt = len(aBas)
    sLVinv = linalg.inv(sLV)
    # Convert aBas to DIRECT COORDINATES
    for iAt in range(nAt):
        aBas[iAt] = matmul(aBas[iAt],sLVinv) # Put positions into
        # "direct" coordinates This keeps the atomic coordinates inside
        # the first unit cell--- not necessary but aesthetically
        # pleasing.

        while any(aBas[iAt] >= (1.0-eps)) or any(aBas[iAt] < (0.0-eps)):
            aBas[iAt] = array([i if i < (1.0-eps) else i-1.0 for i in aBas[iAt]])
            aBas[iAt] = array([i if i >= (0.0-eps) else i+1.0 for i in aBas[iAt]])

    return aBas

def _map_enumStr_to_real_space(system_data,structure_data,minkowskiReduce):
    """Maps an enumerated structure back to real space. Returns a
    dictionary containing the real space data.
    :arg system_data: a dictionary containing all the information about the sysytem.
    :arg sturture_data: a dictionary containing the information for this structure.
    :arg minkowskiReduce: logical indicating if basis should be reduced.
    """
    from numpy import matmul, allclose, matrix, array
    
    nD = system_data["nD"]
    n = structure_data["n"]
    # DEFINE the non-zero elements of the HNF matrix
    a = structure_data["HNF"][0][0]
    b = structure_data["HNF"][1][0]
    c = structure_data["HNF"][1][1]
    d = structure_data["HNF"][2][0]
    e = structure_data["HNF"][2][1]
    f = structure_data["HNF"][2][2]

    pBas = system_data["dvecs"]
    S = structure_data["diag"]
    pLV = system_data["plattice"]
    HNF = structure_data["HNF"]
    eps = system_data["eps"]
    L = structure_data["L"]
    # Compute the superlattice vectors
    sLV = matmul(pLV,HNF).tolist()
    # Find the coordinates of the basis atoms
    gIndx = []

    if minkowskiReduce:
        sLV = list(map(list,zip(*_minkowski_reduce_basis(list(map(list,zip(*sLV))),eps))))
        
    # Find each atomic position from the g-space information
    aBas = []
    ic = 0  # Keep track of the number of points mapped so far
    # Loop over the number at sites/parent cell (the d set)it
    for iD in range(0, nD):  
        # For the limits on the loops, see the "interior_points.pdf" write-up
        for z1 in range(a):
            for z2 in range(int((b*z1)/a), int(c+(b*z1)/a)):
                for z3 in range(int(z1*(d-(e*b)/c)/a+(e*z2)/c), int(f+z1*(d-(e*b)/c)/a+(e*z2)/c)):
                    ic +=1
                    # if ic > n: #pragma: no cover
                    #     err("Problem with basis atoms in map_enpStr_to_real_space...")
                    #     exit()
                    # Atomic basis vector in Cartesian coordinates
                    temp = matmul(pLV,[z1,z2,z3]).tolist()
                    temp2 = [temp[i]+pBas[iD][i] for i in range(len(pBas[iD]))]
                    aBas.append(temp2)
                    # Map position into the group
                    greal = matmul(L,[float(z1),float(z2),float(z3)]).tolist() 
                    g = [int(i) for i in greal] # Convert the g-vector from real to integer
                    if not allclose(greal,g,rtol=eps,atol=eps): #pragma: no cover
                        err("map2G didn't work in map_enumStr_to_real_space")
                        exit()
                    # Bring the g-vector back into the first tile
                    g = [g[i]%S[i] for i in range(len(S))] 
                    # gIndx is the index in the configuration string that
                    # tells us which atom type is used at this position

                    gIndx.append((iD)*S[0]*S[1]*S[2]+g[0]*S[1]*S[2]+g[1]*S[2]+g[2])
    if ic != n*nD: #pragma: no cover
        err("ERROR: map_enumStr_to_real_space: Didn't find the correct # of basis atoms")
        exit()

    k = system_data["k"]
    x = []
    for i in range(k):
        x.append(0.0)
        
    labeling = structure_data["labeling"]
    spin = []
    if k % 2 == 0:
        for iAt in range(0, n*nD):
            i = int(labeling[gIndx[iAt]])
            digit = i-k//2 # convert 0..k-1 label to spin variable -k/2..k/2
            x[i] += 1  # Keep track of the concentration of each atom type
            if digit < 0:
                spin.append(digit)
            else:
                spin.append(digit+1) # skip 0 as a spin if k is even
    else:
        for iAt in range(0, n*nD):
            i = int(labeling[gIndx[iAt]])
            spin.append(i-k//2)
            x[i] += 1 # Keep track of the concentration of each atom type

    x = [i/float(n*nD) for i in x]

    space_data = {"sLV": list(map(list,zip(*sLV))), "aBas": aBas, "spin": spin, "gIndx": gIndx, "x": x}

    return space_data

def _minkowski_reduce_basis(IN,eps):
    """Performs a minkowski reduction on the basis atoms."""

    from numpy import allclose, linalg, array, matrix
    from copy import deepcopy

    limit = 10

    if allclose(linalg.det(IN),0.0,rtol=eps,atol=eps):
        raise ValueError("Input basis for 'minkowski_reduce_basis' was not linearly independent")

    OUT = deepcopy(IN)

    # Keep applying the greedy algorithm until the vectors come out already sorted
    for it in range(1, limit +1):
        # Sort the three vectors into ascending order
        temp = deepcopy(OUT)
        norms = linalg.norm(temp,axis=1).tolist()
        tt = list(range(3))
        tt.reverse()
        for i in tt:
            idx = norms.index(max(norms))
            temp[i] = OUT[idx]
            norms[idx] = 0

        OUT = deepcopy(temp) # Copy the sorted vectors back to OUT
        (OUT[0], OUT[1], OUT[2]) = _reduce_C_in_ABC(OUT[0],OUT[1],OUT[2],eps)
        if linalg.norm(OUT[2]) >= (linalg.norm(OUT[1])-eps):
            break

    if not _minkowski_conditions_check(OUT,eps): #pragma: no cover
        err("ERROR in minkowski_reduce_basis: Minkowski conditions not met."
            "Number of iterations: {}".format(str(limit)))
        exit()
    
    # we want to make sure that the det is positive.
    # NOTE: This *destroys* the mathematical picture of a "greedy reduced basis" (Minkowski), but
    #       from a physical point of view we don't care ;-)
    #       Either way, the basis is as orthogonal as possible.
    if linalg.det(OUT) < 0:
        temp[0] = OUT[1]
        OUT[1] = OUT[2]
        OUT[2] = temp[0]

    return OUT

def _read_enum_out(args):
    """Reads the struct_enum.out file and builds a dictionary with the needed
    information to construct a POSCAR.
    :arg args: The makeStr.py input arguments
    """

    from numpy import transpose
    # which structures are wanted
    if args["structures"] == None:
        with open(args["input"],"r") as f:
            for count, l in enumerate(f):
                pass
        structures = list(range(1,count-13))
    else:
        structures = args["structures"]
    # open the enum.out style file.
    structf = open(args["input"],"r")

    # we'll build a dictionary of the system data and a list of
    # dictionaries for the structures that are wanted.
    structure_data = []
    system = {}
    system["plattice"] = []
    system["dvecs"] = []
    line_count = 1
    system["nD"] = 0
    adjust = 0
    for line in structf:
        temp = line.rstrip()
        if not temp.startswith("#") and "#" not in temp.split()[0]:
            if line_count == 1:
                system["title"] = temp
            if line_count == 2:
                system["bulksurf"] = temp
            if line_count in [3,4,5]:
                vec = [float(v) for v in temp.split() if RepresentsFloat(v)]
                system["plattice"].append(vec)
            if line_count == 6:
                system["nD"] = int(temp.rstrip().split()[0])
            if system["nD"] != 0 and line_count in range(7,7+system["nD"]):
                vec = [float(v) for v in temp.split() if RepresentsFloat(v)]
                system["dvecs"].append(vec)
            if line_count == 7+system["nD"]:
                system["k"] = int(temp.split('-')[0].strip())
            if line_count == 9 + system["nD"]:
                system["eps"] = float(temp.strip().split()[0])
            if line_count == 11 + system["nD"]:
                if "T" == temp.strip().split()[0]:
                    adjust = system["nD"] + system["k"] + 1
                else:
                    adjust = system["nD"]
            if line_count - (14 + adjust) in structures and adjust!=0:
                data = temp.split()
                this_struct = {}
                this_struct["strN"] = int(data[0])
                this_struct["hnfN"] = int(data[1])
                this_struct["hnf_degen"] = int(data[2])
                this_struct["lab_degen"] = int(data[3])
                this_struct["tot_degen"] = int(data[4])
                this_struct["sizeN"] = int(data[5])
                this_struct["n"] = int(data[6])
                this_struct["pgOps"] = int(data[7])
                this_struct["diag"] = [int(data[8]),int(data[9]),int(data[10])]
                this_struct["HNF"] = [[int(data[11]),0,0],[int(data[12]),int(data[13]),0],
                                      [int(data[14]),int(data[15]),int(data[16])]]
                this_struct["L"] = [[int(data[17]),int(data[18]),int(data[19])],
                                    [int(data[20]),int(data[21]),int(data[22])],
                                      [int(data[23]),int(data[24]),int(data[25])]]
                this_struct["labeling"] = data[26]
                if len(data) == 28:
                    this_struct["directions"] = data[27]
                else:
                    this_struct["directions"] = '0'*len(this_struct["labeling"])
                structure_data.append(this_struct)
            line_count += 1
            
    system["plattice"] = transpose(system["plattice"])
    
    return (system, structure_data)

def _write_config(system_data,space_data,structure_data,args,mapping=None):
    """Writes a MTP config style file for the input structure and system
    data.
    :arg system_data: a dictionary of the system_data
    :arg space_data: a dictionary containing the spacial data
    :arg structure_data: a dictionary of the data for this structure
    :arg args: Dictionary of user supplied input.
    :arg mapping: A dictionary of the species mappings if to be used.
    """

    from numpy import array, dot, transpose
    from random import uniform
    from numpy.random import randint
    
    filename = args["outfile"]

    # Get the labeling, group index, structure number and arrow labels
    # from the input data structure.
    labeling = structure_data["labeling"]            
    gIndx = space_data["gIndx"]
    arrows = structure_data["directions"]
    struct_n = structure_data["strN"]

    # The arrow basis.
    arrow_directions = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]
    directions = []

    # Construct the concentrations of the atoms from the labeling by
    # counting the number of each type of atom present in the
    # labeling.
    concs = []
    for i in range(system_data["k"]):
        this_conc = 0
        for atom in range(structure_data["n"]*system_data["nD"]):
            if labeling[gIndx[atom]] == str(i):
                this_conc += 1
        concs.append(this_conc)
    def_title = "{}{}\n".format(str(system_data["title"]),str(structure_data["strN"]))

    # Get the lattice parameter for the atomic species provided by the
    # user.
    if mapping is not None and len(args["species"]) != len(mapping):
        species = [args["species"][i] for i in mapping.values()]
    else:
        species = args["species"]
    
    lattice_parameter, title = _get_lattice_parameter(species,concs,
                                                     system_data["plattice"],system_data["nD"],
                                                      def_title,remove_zeros=True)

    # Find out the directions for each arrow.
    for arrow in arrows:
        directions.append(array(arrow_directions[int(arrow)]))

    sLV = list(lattice_parameter*array(space_data["sLV"]))

    # Start writing the data to the file.
    with open(filename,"a+") as poscar:
        # First write the title and the lattice parameter.
        poscar.write("BEGIN_CFG\n Size\n")
        poscar.write("    {}\n".format(structure_data["n"]*system_data["nD"]))
        poscar.write(" SuperCell\n")
        # Then write out the lattice vectors.
        for i in range(3):
            poscar.write("   {}\n".format("      ".join(
                ["{0: .6f}".format(j) for j in sLV[i]])))
        
        poscar.write("  ")

        poscar.write(" AtomData:  id type       cartes_x      cartes_y      cartes_z\n")
        for iAt in range(structure_data["n"]*system_data["nD"]):
            for ilab in range(system_data["k"]):
                rattle = uniform(-args["rattle"],args["rattle"])
                displace = directions[iAt]*args["displace"]*lattice_parameter
                # If the displacement is non zero and we're `rattling`
                # the system then we need to modify the displacement
                # by the amount being rattled.
                displace += displace*rattle
                if labeling[gIndx[iAt]] == str(ilab):
                    # The final atomic position is the position from
                    # the basis plus the total displacement.
                    out_array = list(array(dot(transpose(sLV),space_data["aBas"][iAt])) + displace)
                    if mapping is None:
                        out_lab = ilab
                    else:
                        out_lab = mapping[ilab]
                    poscar.write("             {0}    {1}       {2}\n".format(iAt+1, out_lab, "  ".join(["{0: .8f}".format(i) for i in out_array])))
        poscar.write(" Feature   conf_id  {}\n".format(title.split()[0]))
        poscar.write("END_CFG\n\n")

def _write_POSCAR(system_data,space_data,structure_data,args):
    """Writes a vasp POSCAR style file for the input structure and system
    data.
    :arg system_data: a dictionary of the system_data
    :arg space_data: a dictionary containing the spacial data
    :arg structure_data: a dictionary of the data for this structure
    :arg args: Dictionary of user supplied input.
    """

    from numpy import array
    from random import uniform

    # Get the output file name.
    if "{}" in args["outfile"]:
        filename = args["outfile"].format(str(structure_data["strN"]))
    else:
        filename = args["outfile"] + ".{}".format(str(structure_data["strN"]))

    # Get the labeling, group index, structure number and arrow labels
    # from the input data structure.
    labeling = structure_data["labeling"]            
    gIndx = space_data["gIndx"]
    arrows = structure_data["directions"]
    struct_n = structure_data["strN"]

    # The arrow basis.
    arrow_directions = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]
    directions = []

    # Construct the concentrations of the atoms from the labeling by
    # counting the number of each type of atom present in the
    # labeling.
    concs = []
    for i in range(system_data["k"]):
        this_conc = 0
        for atom in range(structure_data["n"]*system_data["nD"]):
            if labeling[gIndx[atom]] == str(i):
                this_conc += 1
        concs.append(this_conc)
    def_title = "{} str #: {}\n".format(str(system_data["title"]),str(structure_data["strN"]))

    # Get the lattice parameter for the atomic species provided by the
    # user.
    lattice_parameter, title = _get_lattice_parameter(args["species"],concs,
                                                      system_data["plattice"],system_data["nD"],
                                                      def_title,remove_zeros=args["remove_zeros"])

    # Find out the directions for each arrow.
    for arrow in arrows:
        directions.append(array(arrow_directions[int(arrow)]))

    sLV = space_data["sLV"]

    # Start writing the data to the file.
    with open(filename,"w+") as poscar:
        # First write the title and the lattice parameter.
        poscar.write(title)
        poscar.write("{0:.2f}\n".format(lattice_parameter))
        # Then write out the lattice vectors.
        for i in range(3):
            poscar.write(" {}\n".format(" ".join(
                ["{0: .8f}".format(j) for j in sLV[i]])))
        
        poscar.write("  ")

        # Write the concentrations to the output file. If the species
        # strings were passed in by the user and the user requests
        # there be no zeros in the concentration string then we should
        # remove them from the file. Otherwise we default to leaving
        # them in.
        if (args["remove_zeros"] and args["species"] is not None):
            for ic in concs:
                if ic != 0:
                    poscar.write("{}   ".format(str(ic)))
        else:
            for ic in concs:
                poscar.write("{}   ".format(str(ic)))                    
        
        poscar.write("\n")
        poscar.write("D\n")
        # Now write out the atomic positions to the file.
        for ilab in range(system_data["k"]):
            for iAt in range(structure_data["n"]*system_data["nD"]):
                rattle = uniform(-args["rattle"],args["rattle"])
                displace = directions[iAt]*args["displace"]*lattice_parameter
                # If the displacement is non zero and we're `rattling`
                # the system then we need to modify the displacement
                # by the amount being rattled.
                displace += displace*rattle
                if labeling[gIndx[iAt]] == str(ilab):
                    # The final atomic position is the position from
                    # the basis plus the total displacement.
                    out_array = array(space_data["aBas"][iAt]) + displace
                    poscar.write(" {}\n".format(
                        "  ".join(["{0: .8f}".format(i) for i in out_array.tolist()])))
        
def _make_structures(args):
    """Makes a VASP POSCAR file for the desired structures."""

    (system, structure_data) = _read_enum_out(args)

    # for each structure write the vasp POSCAR
    for structure in structure_data:
        # space_data is a dictionary containing the spacial data for
        # the structure
        space_data = _map_enumStr_to_real_space(system,structure,args["mink"])

        space_data["aBas"] = _cartesian2direct(space_data["sLV"],
                                              space_data["aBas"],system["eps"])

        if args["config"]=="f":
            _write_POSCAR(system,space_data,structure,args)
        elif args["config"]=="t":
            _write_config(system,space_data,structure,args,args["mapping"])
        
def examples():
    """Print some examples on how to use this python version of the code."""
    script = "makeStr: Makes a vasp style POSCAR for the desired system."
    explain = ("For all the examples bellow it is assumed you have already "
               "run the enumeration code and produced an struct_enum.out style file.")
    contents = [("Make a single POSCAR file",
            "To make a POSCAR file for a specific structure listed in the "
            "`struct_enum.out` style file you will need to identify the structure \n number "
            "(the first number of each row in the file) for the structure you want "
            ". For example to make a POSCAR for structure number 10 \n from an `struct_enum.out` "
            "file.","makeStr.py 10 \n"),
           ("Make multilpe POSCARS at once",
            "To make multiple POSCARS for a range of values in the `struct_enum.out` style "
            "file simply list the starting and ending structure numbers \n of the range. "
            "To make POSCARS for every structure in the output file use the word `all`.",
            "makeStr.py 10 20 \n  makeStr.py all \n"),
            ("Find the lattice parameter for the system",
             "To have makeStr.py predict the lattice parameter for the system using "
             "Vegard's Law use the -species option followed by a space \n separated list "
             "of the elements in the system.","makeStr.py 10 -species Al Cu \n"),
            ("Include displacements in POSCAR",
             "If `arrows` (displacement directions) were included in the enumeration "
             "then it is possible to displace them off the lattice points \n when making the "
             "POSCARS using the -displace option followed by the displacement amount "
             "expressed in terms of the lattice parameter. \n In other words if `a` is "
             "the lattice parameter and the atoms were to be displaced by `a/2` then "
             "the command would be:","makeStr.py 10 -displace 0.5 \n"),
            ("Make displacements have different lengths in POSCAR",
             "If `arrows` were included in the model and the `-displace` flag is being "
             "used it is possible to 'rattle' the displacements so that \n they are not all "
             "the same length. Using the `-rattle` option applies a random distribution "
             "to the displacements with the larges change \n in the displacements specified "
             "by the user as a fraction of the displacement given. So if a displacement of "
             "0.5 was given and the \n displacements were to be randomized by 1/4 of that total "
             "displacement the the command would be:",
             "makeStr.py 10 -displace 0.5 -rattle 0.25")]
    required = ("REQUIRED: A `struct_enum.out` file.")
    output = ("RETURNS: A vasp style POSCAR labeled vasp.* where the `*` is replaced "
              "with the structure number for the `struct_enum.out` file.")
    details = ("")
    outputfmt = ("")

    example(script, explain, contents, required, output, outputfmt, details)

script_options = {
    "structures": dict(nargs="+",
                        help=("The desired structure numbers from the struct_enum.out file. This "
                              "can be either a single value or a desired range indicated by "
                              "the starting and stopping structure numbers.")),
    "-displace": dict(default=0.0, type=float,
                        help=("The displacement amount for the arrows in units of the lattice "
                               "parameter. Default is 0.")),
    "-input": dict(default="struct_enum.out",type=str,
                        help=("Override the default 'struct_enum.out' file name.")),
    "-mink": dict(default="t", choices=["t","f"],
                        help=("Sets flag to perform minkowski reduction of the basis (T/F)."
                              " Default is True.")),
    "-species": dict(default=None, nargs="+",type=str,
                        help=("Specify the atomic species present in the system.")),
    "-species_mapping": dict(default=None, nargs="+",type=int,
                        help=("Specify the atomic species numbers for the system. This option "
                              "is only used for making MTP config files.")),
    "-outfile": dict(default="vasp.{}",type=str,
                        help=("Override the default output file names: 'vasp.{structure#}'" 
                              "for the structures.")),
    "-rattle": dict(default=0.0, type=float,
                        help=("Randomizes the positions of the atoms in the POSCAR by no "
                              "more than the fraction of the displacement provided.")),
    "-config" : dict(default="f",choices=["t","f"],
                   help=("make an MTP config file instead of a VASP POSCAR. If the "
                         "MTP config file already exists it will be appended to, not "
                         "overwritten.")),
    "-remove_zeros" : dict(default="f",choices=["t","f"],
                   help=("Remove the zeros from the concentrations string in the 'POSCAR'."))    
}
"""dict: default command-line arguments and their
    :meth:`argparse.ArgumentParser.add_argument` keyword arguments.
"""

def _parser_options():
    """Parses the options and arguments from the command line."""
    #We have two options: get some of the details from the config file,
    import argparse

    pdescr = "POSCAR contstruction."
    parser = argparse.ArgumentParser(parents=[bparser], description=pdescr)
    for arg, options in script_options.items():
        parser.add_argument(arg, **options)
        
    args = exhandler(examples, parser)
    if args is None:
        return

    return args #pragma: no cover

def run(args):
    """Generates the vasp output file for the desired structure.
    """

    if args == None:
        exit()
    if args["structures"] != None :
        if not RepresentsInt(args["structures"][0]) and args["structures"][0].lower() == "all":
            args["structures"] = None
        elif len(args["structures"])  == 1 and RepresentsInt(args["structures"][0]):
            args["structures"] = [int(args["structures"][0])]
        elif len(args["structures"]) == 2:
            args["structures"] = list(range(int(args["structures"][0]),
                                            int(args["structures"][1])+1))
        else:
            raise ValueError("Please enter a single structure number, two structures that "
                             "indicate the first and last structure to be used in the input "
                             "file, or all. The values {} don't match this "
                             "format.".format(args["structures"]))
    else:
        raise ValueError("Please enter a single structure number, two structures that "
                         "indicate the first and last structure to be used in the input "
                         "file, or all. The values {} don't match this "
                         "format.".format(args["structures"]))
    
    if args["species_mapping"] is not None:
        args["mapping"] = {}
        for i in range(len(args["species_mapping"])):
            args["mapping"][i] = args["species_mapping"][i]
    else:
        args["mapping"] = None
        
    if args["remove_zeros"] == "t":
        args["remove_zeros"] = True
    else:
        args["remove_zeros"] = False

    _make_structures(args)
        
if __name__ == '__main__':
    run(_parser_options())
