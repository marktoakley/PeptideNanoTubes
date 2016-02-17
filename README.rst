Description
===========
Software to analyse hydrogen bond patterns in cyclic peptide nanotubes.

These programs count parallel and antiparallel hydrogen bonds in the structures
of cyclic peptide nanotubes. The hydrogen bonds are defined using the
`DSSP<http://dx.doi.org/10.1002/bip.360221211>`
method. The input structures can come from PDB files or PATHSAMPLE databases

Installation
============
To build this software you will need a fortran compiler (I use gfortran).
Build with:
 $ cd src
 $ make

Usage
=====

TubeHbond
---------
Analyses the hydrogen bonding pattern in a PDB file.
Usage:
 TubeHbond <string file>  <int residues> <int rings>
An example PDB file is in the examples directory.

PSHbond
-------
Analyses the hydrogen bonding pattern in a PATHSAMPLE database.
Usage:
 PSHBond <int residues> <int rings> <int res_size>


Reference
=========
If you use this software, please cite:

Mark T. Oakley and Roy L. Johnston, J. Chem. Theory Comput., 2014, 10, 1810-1816.
http://dx.doi.org/10.1021/ct500004k
