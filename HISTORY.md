# Revision History for `enumlib`

## Revision 1.0.2

- Fixed a bug where if there was more than one point in the
  multilattice the code produced following bug `(d,g)-->(d',g')
  mapping failed in get_rotation_perm_lists`. This was repaired by
  replacing line 611 of derivative_structure_generator.f90 with the
  following code `if (any(dgp==0) .and. (any(dap==0)
  .and. .not. use_arrows)) stop "(d,g)-->(d',g') mapping failed in
  get_rotation_perm_lists"`.

## Revision 1.0.1
- Updated the polya submodule.
- Fixed a bug in tree_class.f90 that was causing a segmentation fault
  in the arrow enumeration. The error was caused because using
  `rotated_arrowing(self%A%layer(d)%perms(perm_i,:))` only works of the
  arrays are of the same length. Replaced this line with
  `self%A%layer(d)%perms(perm_i,rotated_arrowing)` which operates as
  this permutation should.

- Fixed the final color mapping in the arrow enumeration. The location
  of the old coloring in the old coloring was incorect because the
  computation used the sizee of the color map as a metric rather than
  the difference in the number of colors vs the number of arrows. As a
  result self%color_map(coloring(site_i)-size(self%color_map,1),2) -1
  was replaced with
  self%color_map(coloring(site_i)-(self%k-self%narrows)-1,2) -1.

## Revision 1.0.0
- Fixed the new functions so they compile with ifort as well.

- Made the aperms argument of get_rotation_perms_lists optional to fix
  a compiler issue.

- Fixed sort_conts so it places the arrows species at the end of the
  list and moved it into the heapsort interface.

- Implemented a polya.x function in the code that can be compiled
  independently from enum.x. This program prints the polya predictions
  for the cell sizes.

## Revision 0.1.3
- Removed the now redundant files src/itertools.f90 and
  src/classes.f90 and the reference to them in the Makefile.

- Added comments to multiple subroutines in tree_class.f90.

- Updated the arrow enumeration code so that it can handle 2D and 3D enumerations properly.


## Revision 0.1.2

- Added checks on all allocations that were added for the arrowed enumeration and the recursively stabilized enumeration.
- Changed the indexing variable names to be more descriptive.
- Made some of the comments and documentation in the recursively stabilized enum code more helpful.
- Changed write_single_labeling, in labeling_related.f90, to take an optional arrow_labeling argument so that when arrows are present it will write out the arrow labeling as well. Removed write_arrow_labeling from arrow_related.f90.

- Added symlib and polya as submodules of enumlib. Updated the Makefile te reflect the change.

- Fixed the default cutoff for the switch from enum3 to enum4 when site restrictions are present.

- Removed the check for a single element enumeration in enum4. (The code now works even for a single element enumeration without this check being present.)


## Revision 0.1.1

Added the source code for the arrowed enumeration. Most of the changes happened in the recursively_stabilized_enum code in labeling_related.f90 and in the tree_class.f90 model where a number of subroutines had to be subtely altered and a new model called arrow_related.f90 has been added to handle the arrow specific methods.

## Revision 0.1.0

Added the source code for the non-arrowed version of the recursively_stabilized_enum code. This entailed adding a new module `tree_class.f90` and adding methods to labeling_related.f90 (recursively_stabilized_enum and write_single_labeling) and sorting.f90 (sort_concs). The new approach has been implemented in derivative_structure_generator.f90

## Revision 0.0.3

Cleaned up the codes XML documentation so that the unittests don't produce as many warnings. Also made a few minor changes in the code that caused errors when compiling with the strict flag on.

## Revision 0.0.2

Fixed a compiler issue that was being caused by the fact that GroupList and permList were moved from tree_class.f90 to enumeration_types.f90 but were not made public in enumeration types.

## Revision 0.0.1

This commit is to move the "extra" or "auxiliary" files in the src folder to the folder aux_src. These files may be useful for some users but are not necessary for the enumeration code to operate. We also updated the Makefile to reflect these changes, moved some types that aren't in enumeration_types.f90 but were declared inside specific programs so that the unit tests can run more smoothly, and added a README to the src and aux_src files.

## Initial Repository: Revision 0.0.0

The first few commits were to get the repo up to scratch and nice and clean with installation instructions etc. This includes the commit for new revision number 0.0.0. It includes an update of the `*.xml` files defining the unit tests so that they work with the distribution directory's structure.

