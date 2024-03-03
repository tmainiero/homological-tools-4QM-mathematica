# Mathematica Homological Toolbox for the Quantum Mechanic

## License: MIT 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description
This is a set of Mathematica packages built to compute the cohomology of the chain complexes introduced in the paper 'Homological Tools for the Quantum Mechanic' ([arXiv:1901.02011](https://arxiv.org/abs/1901.02011)).

For software that can compute ranks of cohomology components in Octave/Matlab, see: <https://github.com/tmainiero/homological-tools-4QM-octave>.

Unlike the Octave software, these packages are written in a functional style, and also allow one to output explicit representatives of generators of cohomology components. 

An interactive quick start guide/minimal documentation is included in `quick_start.nb`.

Most functions are documented.  (If we are interested in a function `functionOfInterest`, a description can be found by evaluating `??functionOfInterest`.)

The package `BasicStable.wl` provides useful quantum mechanical operations and can be used independently of the remainder of the packages.
The repository (https://github.com/tmainiero/quantum-info) contains only this package.

## How to Download

### Git

`git clone https://github.com/tmainiero/homological-tools-4QM-mathematica.git`


### From the Github web interface
Click that fancy green "Clone or download" button on the top right!
Then "Download ZIP".


### File by File

#### From the Github web interface:
1. Go to the file you want to download and click it to view the contents
2. Locate the "Raw" button (On the top right at the time of writing) and right click.
3. Save as...

Make sure that all files are located in the $Path of your Mathematica directory, if not the same folder!

## Caveat
Unfortunately, this software was written before the GNS and commutant cochain complexes were fully understood.  The techniques used to compute cohomology, present in `CechOpsStable`, are meant to compute Cech *homology* of a *co*-presheaves, rather than Cech cohomology of presheaves---the latter being the technique used in the paper.  One can, however, recover the appropriate cohomology after some degree shifting (and sign corrections) in the form of some wrapper functions present in `StateHomologyStable`.  A solution to this confusion would be a rewrite of the functions in the intermediate package `CechOpsStable`; this might be planned for future versions of this software.

## Wolfram Engine as a Mathematica Alternative
If you do not have Mathematica, there is "free" (as in free beer) alternative: namely the Wolfram Engine, which can be run as a Jupyter Kernel.
See Wolfram research's official installation instructions [here](https://github.com/WolframResearch/WolframLanguageForJupyter).
This [YouTube video](https://youtu.be/p0sXuj9lS2k?si=oHdGmbXaINLnhXmz) (circa Feb 2023) demonstrates some capabilities of this kind of setup, and [this accompanying blog post](https://arundquist.wordpress.com/2023/02/23/mathematica-for-free/) for a Window user's setup.

As far as I know, the Wolfram Engine has the same full computational capabilities as in Mathematica, the only drawback of running it in a Jupyter notebook seems to be the loss of some graphical capabilities of Mathematica notebook front end.
For the purposes of interacting with the code here, nothing serious is lost.

## Why Mathematica?
Mathematica has a severe drawback: it is not free or open-source software.
Among other things, this makes it inaccessible to many.
Yet, the choice was made to work with it for this codebase was for two reasons: 

1. Its tight integration of symbolic and numerical capabilities in a streamlined interface: a somewhat key feature that allows one to explore the kind of ideas in ([arXiv:1901.02011](https://arxiv.org/abs/1901.02011)) in a streamlined fashion.
Alternatives such as SageMath just do not have this tight integration or have more serious drawbacks for implementation that requires both numerical and symbolic capabilities.

2. Mathematica is used widely by high-energy physicists, I'm a high-energy physicist by training and pitched the ideas of ([arXiv:1901.02011](https://arxiv.org/abs/1901.02011)) to high energy physicists first.

In the future this software may be ported to alternatives for pragmatic reasons.

