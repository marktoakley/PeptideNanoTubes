'''
Set up a cyclic peptide nanotube and visualise the structure with matplotlib.
@author: Mark Oakley
'''

import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from tubemaker.nanotube import read_amber_coords, orient_coords, build_tube

if __name__ == "__main__":
    # Locate Amber library
    amber_home=os.environ.get('AMBERHOME')
    lib = amber_home+"/dat/leap/lib/all_amino03.lib"
    # Build nanotube
    res_coords =  read_amber_coords("ALA", lib)
    res_coords = orient_coords(res_coords)
    coords = build_tube(2, 8, res_coords)
    # Visualise structure
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111, projection='3d')
    ax.clear()
    xs = coords[:,0]
    ys = coords[:,1]
    zs = coords[:,2]   
    ax.set_xlabel('x')
    ax.set_ylabel('y') 
    ax.set_zlabel('z')
    ax.set_xlim3d(-10,10)
    ax.set_ylim3d(-10,10)
    ax.set_zlim3d(-10,10)
    ax.scatter(xs,ys,zs,s=10)
    plt.show()
