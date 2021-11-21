#!/usr/bin/env python
import sys
import numpy as np
from rdkit import Chem

def cog(atom_num, mass_arr, coords_arr): # Calculate weighed centre of 
    sum_mass = 0
    sum_coord = np.zeros(3) # Create an array to store sum of each atomic coordinate
    for i in range(atom_num):
        sum_mass = sum_mass + float(mass_arr[i])
        weighted_coord = np.multiply(coords_arr[i], float(mass_arr[i])) # Multiply atomic coordinates by atomic mass
        sum_coord = sum_coord + weighted_coord # Add weighted coordinate to sum_coord
    
    weight_cog = sum_coord / sum_mass # Divide sum of weighted coordinates by sum of atomic mass
    return(weight_cog)

def ele_count(mol, ele): # Count atoms of specified element
    count = 0 # Initialise count
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if (sym == ele): 
        # If symbol matches specified element
            count = count +1
            # Count increases by 1
            
    return(count)
            
    
    