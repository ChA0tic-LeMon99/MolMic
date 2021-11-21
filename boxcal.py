#!/usr/bin/env python
import sys
import numpy as np
from rdkit import Chem
def boxcal(mymol):
    icn = 0
    tlist = [] 
    #declare a list tlist going to use tlist to hold a list of 
    #temporrary variable that will be used in the program for many
    #different purposes
    for atom in mymol.GetAtoms():
        pos = mymol.GetConformer().GetAtomPosition(icn)
        
        tlist.append(pos)
        icn = icn + 1
    coords = np.array(tlist)
    minmax = np.zeros((6,))
    for j in range(3):
        minmax[j] = coords[0][j] # Set min and max as first set of coordinates
        minmax[j + 3] = coords[0][j]
    for i in range(icn): # Iterate through each set of coordinates
        for j in range(3):
            if(coords[i][j] > minmax[j]): minmax[j] = coords[i][j] # First three coordinates are max
            if(coords[i][j] < minmax[j + 3]): minmax[j + 3] = coords[i][j] # Last three are min
    print("Max: {} {} {}".format(minmax[0],minmax[1],minmax[2]))
    print("Min: {} {} {}".format(minmax[3],minmax[4],minmax[5]))
    extent = np.zeros(3) # Initialise extent which stores distance between min and max
    size = 1 # Initialise size which stores molecule size
    for j in range(3):
        extent[j] = minmax[j] - minmax[j+3]
        size = size * extent[j] # Calculate size assuming molecule is cuboid
    size = round(size,3)
    print ("Size: {}\n".format(size))
    
    return (size)
    
def sizediff(querysize,refsize): # Calculate differnec between query and ref size
    if (querysize > refsize):
        diff = querysize - refsize
    else:
        diff = refsize - querysize
    diff = round(diff,3)
    return (diff)
    
    


    
