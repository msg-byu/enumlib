# Revision History for `enumlib`

## Revision 0.0.3

Cleaned up the codes XML documentation so that the unittests don't produce as many warnings. Also made a few minor changes in the code that caused errors when compiling with the strict flag on.

## Revision 0.0.2

Fixed a compiler issue that was being caused by the fact that GroupList and permList were moved from tree_class.f90 to enumeration_types.f90 but were not made public in enumeration types.

## Revision 0.0.1

This commit is to move the "extra" or "auxiliary" files in the src folder to the folder aux_src. These files may be useful for some users but are not necessary for the enumeration code to operate. We also updated the Makefile to reflect these changes, moved some types that aren't in enumeration_types.f90 but were declared inside specific programs so that the unit tests can run more smoothly, and added a README to the src and aux_src files.

## Initial Repository: Revision 0.0.0

The first few commits were to get the repo up to scratch and nice and clean with installation instructions etc. This includes the commit for new revision number 0.0.0. It includes an update of the `*.xml` files defining the unit tests so that they work with the distribution directory's structure.

