##BRIEF DESCRIPTION

## Description of aux_src files

These files contain program's that can be useful for some users and need should be maintained but are not necessary for the main code, enum.x, to run. A description of what each file does follows:

#makestr.f90:
makestr.f90 is used for making POSCARs, it is accessed and used in the src directory by typing
`make makestr.x`

#make2Dplot.f90: This program is used to make a 2D plot of the
enumeration results. It can also be compiled in the src directory by
typing `make 2Dplot.x`. This program relies on splot.f which won't
compile with gfortran. To use this program compile with ifort.

#convert_structure_to_enumformat.f90:
This program is designed to convert a structures.in file to structure.enum.out file. It can be compiled in the src directory using `make convert_structures_to_enumformat.x`

#HNF_counter.f90:

#compare_two_enum_files.f90:
It can be compiled from the src directory using `make compare_enum_files.x`

#cwrapper.f90:
This prgram provides an interface for C++.

#find_structure_in_list.f90:
It can be compiled from the src directory using `make find_structure_in_list.x`
ih
#makePerovStr.f90:
This program reads in a file in struct_enum.out format and writes out a VASP-like structue file. It can be compiled from the src directory using `make makeperovstr.x`

#makeStr2D.f90:
This program reads in a file in struct_enum.out format and writes out a VASP-like structue file for a surface calculation. It can be compiled in the src directory using `make makestr.2d`

#random_lattice_driver.f90:

If you have questions, email or call: gus.hart@gmail.com, 801-422-7444

For in-depth examples, please see the `EXAMPLES` file in the `support` directory.