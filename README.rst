A set of programs to analyse hydeogen bond patterns in cyclic peptide nanotubes
in Amber/GMIN/PATHSAMPLE.

This package includes two components:

TubeHbond
---------
Analyses the hydrogen bonding pattern in a PDB file.
Usage: TubeHbond <string file>  <int residues> <int rings>

PSHbond
-------
Analyses the hydrogen bonding pattern in a PATHSAMPLE database.
Usage: PSHBond <int residues> <int rings> <int res_size>

Hydrogen bonds are defined using the DSSP method:
http://dx.doi.org/10.1002/bip.360221211.
