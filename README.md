##BRIEF DESCRIPTION

This code generates the *derivative superstructures* of a parent
lattice. It works for general lattices, including "multilattices," like
HCP, which have more than one lattice point in the unit cell. The code can enumerate over all concentrations or in a restricted concentration range. 

As of Dec. 2011, both modes (all concentrations vs. restricted
conc. range) can also be used with site restrictions---some atomic
types can be disallowed on some sites in the parent lattice. In other
words, the possible occupation of sites by different atoms does not
have to be the same on each parent lattice site. You could have two
sites in the parent lattice and one could be populated by A or B atoms
and the second site could be populated by B and C atoms, for
example. (But see Ex. 5 in the EXAMPLES file for some caveats.)

A discussion of the algorithm and its applications can be found in the
following publications. If you use this code in a publication, we
would appreciate it if you would cite these papers. Thank you.


[Gus L. W. Hart and Rodney W. Forcade, "Algorithm for generating
derivative structures," Phys. Rev. B 77 224115, (26 June 2008)](http://msg.byu.edu/papers/GLWHart_enumeration.pdf)

[Gus L. W. Hart and Rodney W. Forcade, "Generating derivative
structures from multilattices: Application to hcp alloys,"
Phys. Rev. B 80 014120 (July 2009)](http://msg.byu.edu/papers/multi.pdf)

[Gus L. W. Hart, Lance J. Nelson, and Rodney W. Forcade, "Generating
derivative structures at a fixed concentration," Comp. Mat. Sci. 59
101-107 (March 2012)](http://msg.byu.edu/papers/enum3.pdf)

If you download this code, please send me an email so I can let you
know about bugs/improvements/etc.


##COMPILING THE CODE

To compile the code, checkout the `symlib` repository from github, as
well as this repository.

```
git clone https://github.com/msg-byu/enumlib.git  
git clone https://github.com/msg-byu/symlib.git
```

Both repositories should be at the same level in your folder
hierarchy.  

First, compile `symlib`.   
Go to the `symlib/src` directory:  
```cd symlib/src```

Set an environment variable to identify your fortran compiler:  
[in the bash shell, gfortran compiler]  
```export F90=gfortran```  
[the `Makefile` also recognizes `ifort`]

Then compile using the Makefile:  
`make`

(Alternatively, instead of setting the `F90` environmental variable first, you may just specify the variable during the make: `make F90=gfortran`.)

Next, make the enumeration library  
```cd ../../enumlib/src```  
`make`

Finally, to make a stand-alone executable for enumeration:  
```make enum.x```

If you have trouble compiling, please send me a note and I'll try and
assist you.  gus.hart@gmail.com, 801-422-7444

##HOW TO USE THE CODE:

Example input files are in the `input` ('enumlib/support/input') directory. 

There are also helper programs in the aux_src directory for making VASP POSCARs, pictures of
2-D enumerations, etc.

`make makestr.x`  <-- To make a program for making POSCARs  
`make 2Dplot.x` <-- To make a program that plots 2D enumeration results

If you have questions, email or call: gus.hart@gmail.com, 801-422-7444

