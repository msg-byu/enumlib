# READ THE DOCUMENTATION FOR ENUMLIB

#from io import StringIO, BytesIO, FileIO
import argparse
import io
from numpy import array
import os
from random import uniform
import subprocess as sp
import yaml
from multiprocessing import cpu_count, Process, Manager, Queue

from ase import Atoms
from ase.db import connect
from ase.io import read as read_ase
from ase.build import bulk

from mlippy.makeStr import *
from mlippy.makeStr import _read_enum_out, _map_enumStr_to_real_space, _cartesian2direct, _get_lattice_parameter


def parse_args():
    

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('yaml', nargs='+',
                        help='One or more yaml files ready for enumeration')
    parser.add_argument('--clear', type=bool,
                        help='Clear previous database or not')

    args = parser.parse_args()
    return args


class Enumerate():
    """
    Contains all the parameters necessary to enumerate bulk structures for training an MTP

    TODO: Add capability to handle restrictions on lattice points



    makeStr default arguments:
    default_args = {
      'examples': False,
      'verbose': False,
      'action': 'print',
      'debug': False,
      'structures': None, #['all']
     'displace': 0.0,
      'input': 'struct_enum.out',
      'mink': 't',
      'species': [],
      'species_mapping': [],
    #    'outfile': 'vasp.{}',
      'rattle': 0.0,
      'config': 'f',
      'remove_zeros': 'f'
    }        
    """

    def __init__(self, yaml_path):
        # Setup the enumeration class by reading in a yaml file
        self.read_yaml(yaml_path)

        manager = Manager()
        self.db_return = manager.dict()
        print(self.db_return)

        self.SetupProcess()
        self.clear = False
        
        self.set_db_path()


    def set_db_path(self, db_name='enum.db'):
        self.db_path = os.path.join(self.setup_path, db_name)


    def read_yaml(self, yaml_path):
        """Read the configuration yaml file for the alloy database"""

        with open(yaml_path) as yp:
            yml = yaml.safe_load(yp)

        self.species = list(yml['species'].values())
        print(self.species)
        self.nspecies = len(self.species)
        alloy = ''.join(self.species)


        self.database_path = os.path.join(yml['ROOT'], alloy)
        self.setup_path = os.path.join(self.database_path, 'setup')

        
        self.__dict__.update(yml['enumeration'])
        

    def ClearDatabase(self):
        enum_db_path = os.path.join(self.setup_path, 'enum.db')
        if os.path.isfile(enum_db_path):
            os.remove(enum_db_path)
        
    def SetupProcess(self):
        """
        Using the `ase.build.bulk` module, we calculate the lattice vectors and lattice points 
        that will be used in the enumeration.

        
        """

        self.Process = []
        
        for parent in self.lattice:
            enum_path = os.path.join(self.setup_path, parent)
            struct_enum = os.path.join(enum_path, 'struct_enum.in')
            
            if not os.path.isdir(enum_path):
                os.mkdir(enum_path)

            struct = bulk('X', parent, a=1)
            npoints = len(struct.positions)
            
#            print(parent)

            with open(struct_enum, 'w') as f:
                print('Writing input for:', parent)
                print(parent, file=f)
                print('bulk', file=f)

                for vec in struct.cell: 
                    print(*vec, sep=' ', file=f)

                print(self.nspecies, '-nary case', file=f)
                print('  ', npoints, '  # number of points in lattice', file=f)
                
                for pos in struct.positions:
                    print(*pos, self.restrictions, sep=' ', file=f)
                from math import ceil
                cell_size = [ ceil( dim / npoints ) for dim in self.cell_size ]  # fixes size for multilattices
                print(cell_size)
#                cell_size[1] = int(cell_size[1] / npoints)
                print(*cell_size, sep=' ', file=f)
                print(self.finite_precision, file=f)
                print(self.labelings, file=f)
                
                try:
                    for conc in self.concentrations:
                        print(conc, file=f)
                except TypeError:
                    print('No concentration restrictions applied')
            
            # initizlize the Processes
            #    print(enum_path, "ENUM_PATH")
            #self.ClearDatabase()
            #q = 
            #            q.put(enum_path)
            #            p = Process(target=self._EnumerateStructs, args=(q,))
            #return_dict = self.manager.dict()
            p = Process(target=self._EnumerateStructs, args=(enum_path,self.db_return))

#            q.put(p)
            self.Process.append(p)

            
    def _EnumerateStructs(self, enum_path, db_return):
        """
        Worker process for enumeration and generating the atoms database
        """
        with open(os.path.join(enum_path, 'enum.out'), 'wb') as f:
            # Cleanup before enum.x
            my_env = os.environ.copy()

            sp.run(['rm',
                    'enum.out',
                    'vasp.cfg'],
                   cwd=enum_path,
                   stderr=sp.PIPE,
                   stdout=sp.PIPE,
                   env=my_env
            )
            
            # Run enum.x
            print('Running: enum.x')
            sp.run(['enum.x'], cwd=enum_path, stderr=f, stdout=f, env=my_env)

        # Convert the `struct_enum.out` file into an Atoms database
        #args = self.makeStr
        self.makeStr['species'] = self.species

        structs = self.makeStr['structures']
        if type(structs) is str:
            # test to see if a range has been given
            try:
                bounds = [ int(i) for i in structs.split(' ') ]
                if len(bounds) == 2:
                    bounds[1] = bounds[1] + 1
#                    print(bounds)
                    self.makeStr['structures'] = range(*bounds)
            except:
                list_structs = [int(i) for i in structs.split(',')]
                self.makeStr['structures'] = list_structs

            

        # Build the database and run the makeStr functions
        atoms = make_database(enum_path, self.makeStr)
        self.db_return[enum_path] = atoms
        
        

    def RunEnum(self):
        
        for p in self.Process:
            p.start()

        for p in self.Process:
            p.join()    

        
        db_vals = self.db_return.values()
        with connect(self.db_path) as db:
            for parent in db_vals:
                #print("Adding", parent)
                for struct in parent:
                    db.write(struct)

# Write the structures to an ASE object using the makeStr.py template from writing poscars

def make_database(enum_path, args):
    """Writes a vasp POSCAR style file for the input structure and system
    data.
    :arg system_data: a dictionary of the system_data
    :arg space_data: a dictionary containing the spacial data
    :arg structure_data: a dictionary of the data for this structure
    :arg args: Dictionary of user supplied input.

    """

    os.chdir(enum_path)

    (system, structure_data) = _read_enum_out(args)
    atoms = []

    
    for structure in structure_data:
        # space_data is a dictionary containing the spacial data for
        # the structure

        space_data = _map_enumStr_to_real_space(system,structure,args["mink"])
        space_data["aBas"] = _cartesian2direct(space_data["sLV"],
                                                  space_data["aBas"],
                                                  system["eps"])

        atom = write_ASE(system, space_data, structure, args)
        atoms.append(atom)

    print("ATOMS:", len(atoms))
    #return atoms
    database_path = os.path.join('..', 'enum.db')
    #with connect(database_path, append=True) as db:
    #    for struct in atoms:
    #        #print(struct)
    #        db.write(struct)
    print(database_path)
    
    
    return atoms
#    with connect(database_path) as db:
#        counter = 1
 #       
# #       i = len(db)
#        if i == 0:
#            i = 1#
#
#        print('DB len:', i)
#        
#        while counter < len(atoms):
#            i += 1
#            struct_id = db.reserve(name=i)
#            print(i)
#
#            if struct_id is None:
#                print("ID Reserved. Choosing new ID... ", i)
#                
#            else:
#                print("Writing crystal to database: ", counter, i)
#                db.write(atoms[counter], id=struct_id, name=i)
#                counter += 1
#
#            i += 1


## Define the functions that interpret the `struct_enum.out` file


def _write_ASE_OLD(system_data,space_data,structure_data,args):
    """Writes a string in vasp POSCAR style for the input structure and system
    data and then stores that that data as an ASE Atoms object.
    Code adapted from `Enumlib.makeStr.py`.

    :arg system_data: a dictionary of the system_data
    :arg space_data: a dictionary containing the spacial data
    :arg structure_data: a dictionary of the data for this structure
    :arg args: Dictionary of user supplied input.
    """

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
    def_title = '' #"{} str #: {}".format(str(system_data["title"]),str(structure_data["strN"]))
    # Get the lattice parameter for the atomic species provided by the
    # user.
    lattice_parameter, _title = _get_lattice_parameter(args["species"],concs,
                                                      system_data["plattice"],system_data["nD"],
                                                      def_title,remove_zeros=args["remove_zeros"])
    
    title = ' '.join(_title.strip('\n').split())
    for arrow in arrows:
        directions.append(array(arrow_directions[int(arrow)]))

    sLV = space_data["sLV"]

    # Start writing the data to the file.
    with io.StringIO() as poscar:
        # First write the title and the lattice parameter.
        print(title, file=poscar)
        print(lattice_parameter, file=poscar)

        # Then write out the lattice vectors.
        for i in range(3):
            print(*sLV[i], sep='  ', file=poscar)

        cell = sLV


        # Write the concentrations to the output file. If the species
        # strings were passed in by the user and the user requests
        # there be no zeros in the concentration string then we should
        # remove them from the file. Otherwise we default to leaving
        # them in.
        if (args["remove_zeros"] and args["species"] is not None):
            conc_str = ''
            for ic in concs:
                if ic != 0:
                    poscar.write("{}   ".format(str(ic)))
        
#                print(*concs, sep='  ', file=poscar)

        print("\nD", file=poscar)
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
        
                    #        _poscar = poscar.getvalue()
        _poscar = poscar.getvalue()

    atom = read_ase(io.StringIO(_poscar), format='vasp')
    return atom

def write_ASE(system_data,space_data,structure_data,args):
    """Writes a string in vasp POSCAR style for the input structure and system
    data and then stores that that data as an ASE Atoms object.
    Code adapted from `Enumlib.makeStr.py`.

    :arg system_data: a dictionary of the system_data
    :arg space_data: a dictionary containing the spacial data
    :arg structure_data: a dictionary of the data for this structure
    :arg args: Dictionary of user supplied input.
    """
    # Get the labeling, group index, structure number and arrow labels
    # from the input data structure.

    labeling = structure_data["labeling"]
    gIndx = space_data["gIndx"]
    arrows = structure_data["directions"]
    struct_n = structure_data["strN"]

    # The arrow basis.
    arrow_directions = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]


    concs = []
    for i in range(system_data["k"]):
        this_conc = 0
        for atom in range(structure_data["n"]*system_data["nD"]):
            if labeling[gIndx[atom]] == str(i):
                this_conc += 1
        concs.append(this_conc)
 #   print(system_data["k"])
    def_title = '' #"{} str #: {}".format(str(system_data["title"]),str(structure_data["strN"]))

    # Get the lattice parameter for the atomic species provided by the
    lattice_parameter, _title = _get_lattice_parameter(args["species"],concs,
                                                      system_data["plattice"],system_data["nD"],
                                                      def_title,remove_zeros=args["remove_zeros"])

    cell = array(space_data["sLV"]) * lattice_parameter
    formula = ''.join([ args["species"][i] + str(concs[i]) for i in range(len(args['species'])) ])
    positions = space_data["aBas"]
    
    #initialize the ASE.Atoms object
    crystal = Atoms(formula, pbc=True, positions=positions)
    crystal.set_cell(cell, scale_atoms=True)
    crystal.wrap()

    # Apply a random rattle to the structure
    # will eventually be able to take in all ASE.Atoms.rattle() arguments
    try:
        rattle = float(args["rattle"])
    except:
        if args["rattle"]:    # If user gives `True` use default ASE.Atoms.rattle
            crystal.rattle()

    # Test if user input a cell displacement
    try:
        celldisp = float(args["displace"])
    except:
        celldisp = 0

    return crystal

def parallelize(struct_enum_out, n_threads=1):
    """
    struct_enum_out (str): path to the struct_enum_out file that you want to enumerate
    n_threads (int): Number of threads to run on
    """
    import subprocess as sp

    with open(struct_enum_out, "r") as f:
        n_structs = int(f.readlines()[-1].split()[0])

    structs_per_thread = n_structs // n_threads

    start = 1
    for i in range(n_threads):
        if i+1 == n_threads:
            end = n_structs  # This makes sure to grab the last few structures that weren't assigned with the floor division.
        else:
            end = start + structs_per_thread - 1

    output_file = ''.join(["enum.tmp.", str(i+1)])  # set the output filename
    args = ' '.join(["python", "~/automtp/mlippy/mlippy/",
                     str(start), str(end), "-config", "t", "-outfile", output_file, "&"])
    sp.run(args,
           shell=True,
           cwd="/home/hayden/Documents/msg/examples",
          )
    print("thread", i+1, ":", start, end, "enumerated:", end-start)
    # 


    start = end + 1  # set the next start value to the next structure




def main():
    args = parse_args()

    for alloy in args.yaml:
        print("Enumerating:", alloy)
        enum_alloy = Enumerate(alloy)
#        print(enum_alloy.cell_size, type(enum_alloy.cell_size))
#        enum_alloy.ClearDatabase()
        enum_alloy.RunEnum()



    


if __name__ == "__main__":
    main()
