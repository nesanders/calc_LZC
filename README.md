LZC Relation Metallicity Calculator
===================================

This is a python module to calculate galaxy metallicities given two bands of K-corrected photometry using the LZC relation as calibrated by Sanders et al. 2012 (arXiv:1210:5520, http://arxiv.org/abs/1210.5520).

Bug reports and feedback are welcome, please visit https://github.com/nesanders/calc_LZC

Copyright 2012 Nathan E. Sanders.

Learn more https://www.cfa.harvard.edu/~nsanders/papers/LZC/summary.htm

Installation
============

Download and unpack the soure archive.

    wget https://github.com/nesanders/calc_LZC/zipball/master

Run the install script

    sudo python setup.py install

Usage
=====

Import the module from within a python environment

    from LZC import LZC

Now call the LZC calulator function as

    LZC.cLZC( options )

For help using this function, use

    help('LZC.LZC.cLZC')

Here is an example, calculating the LZC metallicity in the PP04 O3N2 diagnostic using the r band luminosity and r-i color.  This will return the metallicity, the Monte Carlo uncertainty (taking into account photometric errors), the intrinsic scatter in the LZC relation, and a flag integer that equals one if the galaxy is outside the calibrated mu range.

    LZC.cLZC(r=-20,i=-20.1,dr=.1,di=.2,M='r',col='ri',diag='PP04_O3N2')

