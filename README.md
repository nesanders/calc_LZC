LZC Relation Metallicity Calculator
===================================

This is a python module to calculate galaxy metallicities given two bands of K-corrected photometry using the LZC relation as calibrated by Sanders et al. 2012.

Bug reports and feedback are welcome, please visit https://github.com/nesanders/calc_LZC

Copyright 2012 Nathan E. Sanders.

`Learn more <https://www.cfa.harvard.edu/~nsanders/papers/LZC/summary.htm>`_.

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

