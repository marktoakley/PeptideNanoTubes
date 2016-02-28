Description
===========
Software to construct and analyse cyclic peptide nanotubes.

Tubemaker
---------
This program constructs Amber_ input files containing cyclic peptide nanotubes.

TubeHbond and PSHbond
---------------------
These programs count parallel and antiparallel hydrogen bonds in the structures
of cyclic peptide nanotubes. The hydrogen bonds are defined using the
DSSP_ method. The input structures can come from PDB_ files or PATHSAMPLE_ databases

Installation
============
Tubemaker
---------
Tubemaker depends on the numpy python library. It also requires an Amber_ library file (several are available in AmberTools, which can be obtained from the Amber_ website).

TubeHbond and PSHbond
---------------------
To build these programs you will need a fortran compiler (I use gfortran).
Build with::
  cd src
  make

Usage
=====
Tubemaker
---------
Tubemaker constructs an Amber_ inpcrd file containing a cyclic peptide nanotube.
Run Tubemaker with::
  python nanotube.py <#rings> <#residues> <res_name> [arguments]
For example, build an antiparallel tetramer of cyclic octa-alanine with
::
  python nanotube.py 4 8 ALA --anti
The initial coordinates for the peptides are taken from an Amber library file.
By default, the ff03 library from Amber or Ambertools is used (if one of them is installed).
Otherwise, a library file can be specified with the --lib argument.
For a full list of arguments, use::
  python nanotube.py -h

TubeHbond
---------
Analyses the hydrogen bonding pattern in a PDB_ file.
Usage::
  TubeHbond <string file> <int residues> <int rings>
The arguments are the name of a PDB file, the number of residues per ring
and the number of cyclic peptide rings. An example PDB file is in the examples
directory.

PSHbond
-------
Analyses the hydrogen bonding pattern in a PATHSAMPLE_ database.
Usage::
  PSHBond <int residues> <int rings> <int res_size>
The arguments are the number of residues per ring, the number of residues per
ring and the number of atoms per residue. Currently, only cyclic peptides
containing one type of residue are supported (but d- and l- variants of the
same residue are fine).

Reference
=========
If you use this software, please cite:

Mark T. Oakley and Roy L. Johnston, J. Chem. Theory Comput., 2014, 10, 1810-1816.
http://dx.doi.org/10.1021/ct500004k

.. _DSSP: http://dx.doi.org/10.1002/bip.360221211
.. _PDB: http://www.rcsb.org/
.. _PATHSAMPLE: http://www-wales.ch.cam.ac.uk/PATHSAMPLE/
.. _GMIN: http://www-wales.ch.cam.ac.uk/GMIN/
.. _Amber: http://ambermd.org
