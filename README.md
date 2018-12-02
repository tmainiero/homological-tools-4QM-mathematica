# Mathematica Homological Toolbox for the Quantum Mechanic

## License: MIT 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description
This is a set of Mathematica packages built to compute the cohomology of the chain complexes introduced in the paper 'A Homological Toolkit for the Quantum Mechanic.'

Unlike the Octave software, these packages are written in a functional style, and also allow one to output explicit representatives of generators of cohomology components. 

An interactive quick start guide/minimal documentation is included in `quick_start.nb`.

Most functions are documented.  The package `BasicStable.wl` provides useful quantum mechanical operations and can be used independently of the remainder of the packages.

## How to Download

### Git

`git clone https://github.com/tmainiero/homological-tools-4QM-mathematica`

### File by File

#### From the Github web interface:
1. Go to the file you want to download and click it to view the contents
2. Locate the "Raw" button (On the top right at the time of writing) and right click.
3. Save as...

Make sure that all files are located in the $Path of your Mathematica directory, if not the same folder!


## Caveat
Unfortunately, this software was written before the GNS and commutant cochain complexes were fully understood.  The techniques used to compute cohomology, present in `CechOpsStable,` are meant to compute Cech *homology* of a *co*-presheaves, rather than Cech cohomology of presheaves---the latter being the technique used in the paper.  One can, however, recover the appropriate cohomology after some degree shifting in the form of some wrapper functions present in `StateHomologyStable`.  A solution to this confusion would be a rewrite of the functions in the intermediate package `CechOpsStable`; this might be planned for future versions of this software.
