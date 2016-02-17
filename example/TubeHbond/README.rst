TubeHbond
+++++++++

Usage:
 $ TubeHbond <string file>  <int residues> <int rings>

Example: TubeHBond tetramer.pdb 8 4
Counts the hydrogen bonds in the included pdb file (a tetramer of a cyclic
octapeptide). This is a perfect antiparallel nanotube, so you should get 24
antiparallel hydrogen bonds and no parallel ones.
