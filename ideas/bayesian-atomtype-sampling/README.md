# Bayesian atom type sampling example

This is a simple example of how Bayesian atom type sampling using reversible-jump Markov chain Monte Carlo (RJMCMC) [1] might work.

## Manifest

```
atomtype-sampling-example.py - example illustrating use of RJMCMC to sample over SMARTS-specified atom types
datasets/ - small molecule datasets for hydration free energy fitting example
atomtypes/ - input files for atom type sampling
```

## Documentation

### How it works

Atom types are specified by SMARTS matches with corresponding parameter names.

First, we start with a number of initial "base types", specified in `atomtypes/basetypes.smarts`:
```
% atom types
[#1]    hydrogen
[#6]    carbon
[#7]    nitrogen
[#8]    oxygen
[#9]    fluorine
[#15]   phosphorous
[#16]   sulfur
[#17]   chlorine
[#35]   bromine
[#53]   iodine
```
Note that lines beginning with `%` are comment lines.

Atom type creation moves attempt to split off a new atom type from a parent atom type by combining (via an "and" operator, `&`) the parent atom type with a "decorator".
The decorators are listed in `atomtypes/decorators.smarts`:
```
% bond order
$([*]=[*])     double-bonded
$([*]#[*])     triple-bonded
$([*]:[*])     aromatic-bonded
% bonded to atoms
$(*~[#1])      hydrogen-adjacent
$(*~[#6])      carbon-adjacent
$(*~[#7])      nitrogen-adjacent
$(*~[#8])      oxygen-adjacent
$(*~[#9])      fluorine-adjacent
$(*~[#15])     phosphorous-adjacent
$(*~[#16])     sulfur-adjacent
$(*~[#17])     chlorine-adjacent
$(*~[#35])     bromine-adjacent
$(*~[#53])     iodine-adjacent
% degree
D1             degree-1
D2             degree-2
D3             degree-3
D4             degree-4
D5             degree-5
D6             degree-6
% valence
v1             valence-1
v2             valence-2
v3             valence-3
v4             valence-4
v5             valence-5
v6             valence-6
% total-h-count
H1             total-h-count-1
H2             total-h-count-2
H3             total-h-count-3
% aromatic/aliphatic
a              atomatic
A              aliphatic
```
Each decorator has a corresponding string token (no spaces allowed!) that is used to create human-readable versions of the corresponding atom types.

For example, we may find the atom type ```[#6]&H3``` which is `carbon total-h-count-3` for a C atom bonded to three hydrogens.

Newly proposed atom types are added to the end of the list.
MAfter a new atom type is proposed, all molecules are reparameterized using the new set of atom types.
Atom type matching proceeds by trying to see if each SMARTS match can be applied working from top to bottom of the list.
This means the atom type list is hierarchical, with more general types appearing at the top of the list and more specific subtypes appearing at the bottom.

If a proposed type matches zero atoms, the RJMCMC move is rejected.

Currently, the acceptance criteria does not include the full Metropolis-Hastings acceptance criteria that would include the reverse probability.  This needs to be added in.

## Prerequisites

* Python
[Anaconda](https://store.continuum.io/cshop/anaconda/) is suggested.

* OpenMM
Suggested installation route is via [conda](http://conda.pydata.org):
```
conda install openmm
```


## Usage

Example:
```
python atomtype-sampling-example.py --basetypes=atomtypes/basetypes.smarts --decorators=atomtypes/decorators.smarts --molecules=datasets/mobley-504-molecules.sdf --iterations 1000
```

Initially, the base atom types are added to the pool of current atom types, and the number of atoms and molecules matched by each atom type are shown:
```
INDEX        ATOMS  MOLECULES                                          TYPE NAME                                           SMARTS
    1 :       4371        483 |                                         hydrogen                                             [#1]
    2 :       2717        488 |                                           carbon                                             [#6]
    3 :        117        102 |                                         nitrogen                                             [#7]
    4 :        338        216 |                                           oxygen                                             [#8]
    5 :         69         24 |                                         fluorine                                             [#9]
    6 :          2          2 |                                      phosphorous                                            [#15]
    7 :         20         18 |                                           sulfur                                            [#16]
    8 :        116         62 |                                         chlorine                                            [#17]
    9 :         28         23 |                                          bromine                                            [#35]
   10 :         12         11 |                                           iodine                                            [#53]
TOTAL         7790        490
```
After many iterations, the pool of current atom types will have diverged, with some children having been added to the set.  (Base atom types can never be deleted.)
```
Iteration 999 / 1000: True
INDEX        ATOMS  MOLECULES                                          TYPE NAME                                           SMARTS
    1 :       4366        483 |                                         hydrogen                                             [#1]
    2 :          5          5 |                         hydrogen sulfur-adjacent                                  [#1&$(*~[#16])]
    3 :       2602        474 |                                           carbon                                             [#6]
    4 :         22         18 |                         carbon fluorine-adjacent                                   [#6&$(*~[#9])]
    5 :          7          6 |       carbon fluorine-adjacent hydrogen-adjacent                         [#6&$(*~[#9])&$(*~[#1])]
    6 :         25         23 |                          carbon bromine-adjacent                                  [#6&$(*~[#35])]
    7 :         61         33 |         carbon total-h-count-1 nitrogen-adjacent                                [#6&H1&$(*~[#7])]
    8 :        105         92 |                                         nitrogen                                             [#7]
    9 :         12         12 |                           nitrogen triple-bonded                                  [#7&$([*]#[*])]
   10 :        338        216 |                                           oxygen                                             [#8]
   11 :         69         24 |                                         fluorine                                             [#9]
   12 :          2          2 |                                      phosphorous                                            [#15]
   13 :         16         14 |                                           sulfur                                            [#16]
   14 :          4          4 |                                 sulfur valence-6                                         [#16&v6]
   15 :        116         62 |                                         chlorine                                            [#17]
   16 :         28         23 |                                          bromine                                            [#35]
   17 :         12         11 |                                           iodine                                            [#53]
TOTAL         7790        490
```

## References

[1] Green PJ. Reversible jump Markov chain Monte Carlo computation and Bayesian model determination. Biometrika 82:711, 1995.
http://dx.doi.org/10.1093/biomet/82.4.711

