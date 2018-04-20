# Notes about future developments for enumlib

April 19 2019

## Missing use cases
Some users have attempted to use enumlib to do
enumerations that were not part of the use cases orginally
anticipated. The resulted in the third and fourth enumeration-related
papers, but users continue to come up with scientifically sound,
reasonable use cases that are difficult or impractical for the current
enumeration code in the `enumlib` repo.

One class of missing use cases includes problems where one group of
sites in the unit cell have no configurational degree of freedom, that
is, there is only once kind of atom allowed on those
sites. Combinatorically speaking this is trival (like multiplying by
one) but with the current data structures and hash tables this leads
to significant overhead, even making some problems completely
impractical that should be possible.

A closely related use case consists of more than one group of sites
that are mutually exclusive. For example, imagine decorating a
two-lattice that can take atoms of type A,B on the first site and
types C,D on the second site. This case doesn't have the full
complexity of a quaternary case---it's more like to "coupled" binary
cases---but the current algorithms are drastically slowed down by the
"quarternary-ness" of the problem.

In the preceding example, the combinatorial possibilities of site 1
have no bearing on the possibilities of site 2, so the two problems
could be done separately and then the final answer would by a kind of
"outer product" of the two results. Instead of a 4^k type problem, we
have two 2^(k/2) type problems. 

## Review of the development up to this point

The original enumeration algorithm (papers I and II) addressed the
problem of enumerating colorings of a lattice where the full
concentration range of colors was allowed and every color was allowed
on every site. This was not efficient for cases where the
concentrations were fixed (or limited to a small range) or where some
colors were to allowed on some sites. In principle, one could do the
"full" enumeration problem and discard those configurations that
violated site restrictions or that were outside the desired
concentrations. But with a little thinking it's clear that this is
woefully inefficient in many cases.

In enumeration paper III (CMS, 2012), a new hashing algorithm was
developed that allowed the enumeration to be limited to a single
concentration, but it utilized the same concepts as in the original
algorithm. At this stage site restrictions were also introduced, but
in a simple way---as each configuration was generated, it was checked
to see if it violated site restrictions. If it did, that configuration
(and all those below it in the "tree") was skipped.

Enumeration paper IV (CMS, 2017) essentially used the same concept as III for skipping configurations that violated site restrictions---it built the configurations one color at a time (instead of building the entire configuration) skipping "partial" colorings that were
equivalent to those aleady seen. That allowed "early exits" from
branches of the tree search, skipping duplicates of earlier
configurations in the tree. 

## Ideas for a refactor
* It's impractical in many problems to store
(even just the hashes of) all surviving configurations at the end of
the enumerations. Better to print them out in batches as the
enumeration proceeds. So we need data structures that allow this.
* It seems that if we have disjoint sets of types, we can enumerate the different sites separately, for each set, then do a kind of "outer product" to generate the full list of configurations. In general, doing several smaller problems, separately and then combining should be much, much, more effecient because of the combinatorial explosion that happens with additional colors.
* It would be nice to do a 2D problem first, maybe a few toy problems, to build intuition. Nice undergraduate problem.
* There might be some good ideas in this paper DOI 10.1186/s13321-016-0129-3, Supercell program: a combinatorial structure-generation approach for the local-level modeling of atomic substitutions and partial occupancies in crystals Kirill Okhotnikov, Thibault Charpentier, and Sylvian Cadars
* 


