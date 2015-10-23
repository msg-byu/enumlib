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

Example input files are in the `input` (enumlib/support/input) directory. 

There are also helper programs for making VASP POSCARs, pictures of
2-D enumerations, etc.

`make makestr.x`  <-- To make a program for making POSCARs  
`make 2Dplot.x` <-- To make a program that plots 2D enumeration results

If you have questions, email or call: gus.hart@gmail.com, 801-422-7444

For in-depth examples, please see the `EXAMPLES` file in the `support` directory.