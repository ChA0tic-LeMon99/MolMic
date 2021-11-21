#!/usr/bin/env python
import sys
import numpy as np
from rdkit import Chem

def coords(mol): # Obatin a 2D array of atomic coordinates in a molecule  
    icn = 0
    tlist = [] 
    #declare a list tlist going to use tlist to hold a list of 
    #temporrary variable that will be used in the program for many
    #different purposes
    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(icn)        
        tlist.append(pos)
        icn = icn + 1
    coords = np.array(tlist)     
    
    return (coords)
    
def name(mol): # Obtain name of molecule
    name = mol.GetProp("_Name")
    
    return(name)

def ele(mol): # Take a molecule as argument to obtain a list of elements 
    ele = []
    for atom in mol.GetAtoms():
        ele.append(atom.GetSymbol())
     
    return (ele)
    
def mass(mol): # Take a molecule as argument to obtain an array of element masses 
    tlist = []
    for atom in mol.GetAtoms():
        tlist.append(atom.GetMass())
    mass = np.array(tlist)
     
    return (mass)
    
def smile(mol): # Obtain molecule description in smile format 
    Chem.MolToSmiles(mol)


    

        
    

    