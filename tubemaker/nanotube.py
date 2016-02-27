'''
tubemaker
This program constructs cyclic peptide nanotubes
@author: Mark Oakley
'''

import numpy as np
import os
from math import cos,sin,radians,pi

def read_amber_coords(residue, amber_lib):
    '''Read the coordinates of a peptide residue from an Amber library file.
    Parameters:
    residue- The three-letter amino acid code of a residue (e.g. ALA for alanine)
    amber_lib- The location of an Amber library file'''
    header = "!entry."+residue+".unit.positions"
    coords = []
    with open(amber_lib,'r') as amber_file:
        reading_coords=False
        for line in amber_file:
            if header in line: #Start of coordinates record
                reading_coords=True
            elif residue in line: #End of coordinates record
                reading_coords=False
            elif reading_coords:
                coords.append(line.split())
    if coords == []:
        raise ValueError(residue)
    coords_array=np.array(coords,dtype=float)
    return coords_array

def orient_coords(coords_array):
    '''Re-orient the coordinates from an amber library so that:
    The peptide runs along x axis
    The side chain projects along x axis
    The hydrogen bonds point along z axis'''
    #Translate alpha carbon to origin
    ca = coords_array[2,:]
    coords_array = coords_array - ca
    # Rotate coordinates
    theta = radians(35)
    rot = np.array([[cos(theta),0,sin(theta)],
                     [sin(theta),0, -cos(theta)],
                     [         0,         1. ,   0000]])
    coords_array = coords_array.dot(rot)
    return coords_array

def build_ring(num_res, res_coords):
    '''Construct the coordinates of a cyclic peptide ring.
    Parameters:
    num_res: The number of residues in the ring
    res_coords: The coordinates of a peptide residue (see orient_coords)'''
    radius = 3.8*num_res/(2*pi)
    shift = np.array([0,radius,0])
    res_coords = res_coords + shift
    coords = np.copy(res_coords)
    theta = -2.*pi / num_res
    rot = np.array([[cos(theta),sin(theta),0],
                    [-sin(theta),cos(theta),0],
                    [0,0,-1]])
    for i in range(num_res-1):
        res_coords=res_coords.dot(rot)
        coords = np.append(coords,res_coords,axis=0)
    return coords

def build_tube(num_rings, num_res, res_coords, parallel=True, spacing=5.0):
    '''Construct the coordinates of a cyclic peptide ring.
    Parameters:
    num_rings:
    num_res: The number of residues in the ring
    res_coords: The coordinates of a peptide residue (see orient_coords)'''
    ring_coords = build_ring(num_res, res_coords)
    tube_coords = ring_coords.copy()
    shift = np.array([0.,0.,spacing])
    for i in range(num_rings-1):
        ring_coords=ring_coords+shift
        if not parallel:
            ring_coords = ring_coords*np.array([-1.,1.,1.])
            theta = 2*pi/num_res
            rot = np.array([[cos(theta),-sin(theta),0],
                            [sin(theta), cos(theta),0],
                            [         0,  0,1]])
            ring_coords=ring_coords.dot(rot)
        tube_coords = np.append(tube_coords,ring_coords,axis=0)
    return tube_coords
        

def write_amber_coords(coords,res_name,file_name="coords.inpcrd"):
    '''Write an Amber inpcrd file containing the nanotube's coordinates.'''
    with open(file_name,'w') as inpcrd:
        inpcrd.write(res_name+"\n")
        inpcrd.write("%5i\n" % len(coords))
        new_coords = coords.flatten().tolist()
        for i in range(len(new_coords)):
            inpcrd.write("%12.7f" % new_coords[i])
            if (i+1)%6 ==0:
                inpcrd.write("\n")

if __name__ == "__main__":
    amber_home=os.environ.get('AMBERHOME')
    lib = amber_home+"/dat/leap/lib/all_amino03.lib"
    num_res = 8
    num_rings = 4
    res_name="ALA"
    res_coords =  read_amber_coords(res_name, lib)
    res_coords = orient_coords(res_coords)
    #coords = build_ring(num_res, res_coords)
    coords = build_tube(num_rings, num_res, res_coords,parallel = False)
    write_amber_coords(coords,res_name)
