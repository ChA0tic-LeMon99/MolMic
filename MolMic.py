#!/usr/bin/env python
                                   ##==Import modules==##
import sys
import copy
import os.path
import numpy as np
from rdkit import RDConfig
from rdkit import Chem, DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
fName=os.path.join(RDConfig.RDDataDir,'FunctionalGroups.txt')
# Specify a file that contains information on functional groups
from rdkit.Chem import FragmentCatalog

                                     ##==Import custom modules==##
import CoordAnalyse1 as xyz
import sanitise
import boxcal
import molclass as mlc
                                     ##==Funtion definitions==##

def prompt(): # Define a function prompt() that stores output filename
    output_name = input("Enter name for output file\n") 
    # Prompt user for output filename and assigned it to variable output_name
    try:
        f = open("%s.dat" % output_name, "x") 
        # Check if output file already exists        
    except:
        print("Error: Filename already in use")
        prompt() 
        # Use recursion to keep calling prompt() until a valid output filename is entered
    return output_name  
    # Return the value of valid output filename
                                     
                                     ##Object class##
                                     
class mol_class: 
# Create a class to store molecule properties: molecule name, size, element, coordinates, atomis mass
    def __init__(self, name, smile, ele, coords, mass):
        self.name = name
        self.smile = smile
        self.ele = ele
        self.coords = coords
        self.mass = mass
          
class frag_class: 
# Create a class to store molecule fragment properties: molecule name, molecule smiles, molecule fragment number
     def __init__(self, name, smile, frag):
        self.name = name
        self.smile = smile
        self.frag = frag
           
                                     
if(len(sys.argv) != 3) : 
# Check there are 3 command line arguments
    print("Usage: boxcal.py query.sdf reference.sdf")
    sys.exit(-1)
    
queryfile = sys.argv[1] 
reffile = sys.argv[2]
output_name = prompt()

try:
    m = Chem.MolFromMolFile(queryfile) 
    # Read molecule in query file
except:
    print ("Invalid query file")
    sys.exit(-1)
queryname = m.GetProp("_Name") # Get query molecule name
print(queryname)  # Print query molecule name
querymol = mol_class(copy.deepcopy(mlc.name(m)), copy.deepcopy(mlc.smile(m)), copy.deepcopy(mlc.ele(m)), copy.deepcopy(mlc.coords(m)), copy.deepcopy(mlc.mass(m))) 
# Transfer the value for each 'mol_class' property to object 'querymol'
weight_cog = xyz.cog(len(querymol.mass),querymol.mass,querymol.coords)
print ("Weighted cog: ", str(round(weight_cog[0], 3)), ",", str(round(weight_cog[1], 3)), ",", str(round(weight_cog[2], 3)))
querysize = float(boxcal.boxcal(m)) 
# Print minmax, extent, and size of query

f = open("%s.dat" % output_name, "a") 
# Open output file and write query name, cog, and size
f.write("Query: {}\nWeighted cog: {} {} {}\nSize: {}\n".format(querymol.name, round(weight_cog[0], 3), round(weight_cog[1], 3), round(weight_cog[2], 3), querysize))
f.close()
   
 
try:
    refmol = Chem.SDMolSupplier(reffile) 
    # Read reference molecules
except:
    print ("Invalid reference file")
    sys.exit(-1)
    
refcount = 0 
# Count number of reference molecules
refmollist = [] 
# A list for storing refmol objects
sizediff = [] 
# sizediff array to store size difference between ref and query 
sortedsize = [] 
# Array to store sorted size difference 
match = str(refmol[0].GetProp("_Name")) 
# Set size match to first reference molecule
matchdiff = float(boxcal.boxcal(refmol[0]))
for mol in refmol:
    refmol = mol_class(copy.deepcopy(mlc.name(mol)), copy.deepcopy(mlc.smile(mol)), copy.deepcopy(mlc.ele(mol)), copy.deepcopy(mlc.coords(mol)), copy.deepcopy(mlc.mass(mol))) 
    # Transfer the value for each 'mol_class' property to object 'querymol'
    refmollist.append(refmol) # Insert refmol into list of mol_class objects
    weight_cog = xyz.cog(len(refmol.mass),refmol.mass,refmol.coords)
    # Calculate weighted cof for ref molecule
    print(refmol.name)
    print ("Weighted cog: ", str(round(weight_cog[0], 3)), ",", str(round(weight_cog[1], 3)), ",", str(round(weight_cog[2], 3)))
        
    
    refsize = copy.deepcopy(boxcal.boxcal(mol)) 
    diff = copy.deepcopy(boxcal.sizediff(querysize,refsize)) 
    # Calculate query and ref size difference
    sizediff.append(diff) 
    # Store query and ref size difference
    sortedsize = np.copy(sizediff) 
    #Copy sizediff and sort it
    sortedsize.sort()
    f = open("%s.dat" % output_name, "a") 
    # Open output file and write ref name, cog, and size
    f.write("\n{}\nWeighted cog: {} {} {}\nSize: {}\n".format(refmol.name, round(weight_cog[0], 3), round(weight_cog[1], 3), round(weight_cog[2], 3), refsize))
    f.close()
    refcount = refcount + 1
print ("There are {} reference molecules".format(refcount)) # Print number of ref molecules


refcount = 0
matcharr = [] 
# Match array to store top size matches
diffarr = np.zeros(4)
for mol in range(len(refmollist)):
    refmol = refmollist[refcount]
    diff = sizediff[refcount]
    
    for i in range(4): 
    # For the 4 smallest size differences
        if (diff == sortedsize[i]): 
        # If refmol difference matches smallest difference 
            matcharr.insert(i, refmol.name) 
            # Inserts name into match array
            diffarr[i] = diff               
    
    if (refcount == len(sortedsize) - 1):
        for i in range (4):
            print ("Size match {}: {} with difference {}".format(i+1,matcharr[i],diffarr[i]))
            # Print 4 closest size matches
            f = open("%s.dat" % output_name, "a") 
            # Open output file and write 4 closest size matches
            f.write("\nSize match {}: {} with difference {}".format(i+1,matcharr[i],diffarr[i]))
            f.close()
        
    refcount = refcount + 1
                                        
                                    ##Sanitise query file before fragmentation##
                                    
sanitised = str(sanitise.sanitise(queryfile))
                                    
                                    ##Read query and fragment query##

patt = []
fparams = FragmentCatalog.FragCatParams(1,6,fName)
# Use functional groups in fName as parameters
fparams.GetNumFuncGroups()
fcat=FragmentCatalog.FragCatalog(fparams)
fcgen=FragmentCatalog.FragCatGenerator() 
# Fragment query molecule

query = Chem.MolFromMolFile(sanitised)
num_frag = fcgen.AddFragsFromMol(query,fcat)
# Number of fragments
print ("{} fragments are generated".format(num_frag)) 
# Print number of query fragments

fpgen = FragmentCatalog.FragFPGenerator()
fp = fpgen.GetFPForMol(query,fcat)
# Use fragments as fingerprints
andfp = fp&fp
obl = list(andfp.GetOnBits())
# A list of fingerprints
print ("Query: ",query.GetProp("_Name"))
print ("Smile ",Chem.MolToSmiles(query))
# print query smiles
print ("Fragments: ",len(obl),"\n")
#print(fp.GetNumOnBits())

for i in range (num_frag):
    patt.append(fcat.GetEntryDescription(i))
    # Fragment description
    i = i+1


                                ##Read reference and fragment reference##
                            
ref = Chem.SDMolSupplier('lfp.sdf')
# Read reference molecules
reffraglist = []
# list to store ref fragment objects
oblarr = np.zeros(refcount)
# array to store fragment number of each ref molecule
sortoblarr = np.zeros(refcount)
# array for sorted fragment number
i = 0
for mol in ref:    
    fpgen = FragmentCatalog.FragFPGenerator()
    fp2 = fpgen.GetFPForMol(mol,fcat)
    andfp = fp&fp2
    # Fingerprint query and ref
    obl = list(andfp.GetOnBits())
    # A list of matching fragments  
    name = mol.GetProp("_Name")
    # Get ref molecule name
    reffrag = frag_class(name, Chem.MolToSmiles(mol), len(obl))
    # Assign values of frag_class properties to current ref molecule
    oblarr[i] = reffrag.frag
    # Store current fragment number in oblarr array
    
    print ("Molecule: ",reffrag.name, "\n", "Smile: ", reffrag.smile, "\n", "Fragments matched: ", reffrag.frag, "\n")
    # print ref molecule name, smile, and fragment number
    reffraglist.append(reffrag)
    # Add current ref molecule to reffraglist    
    
    i = i +1    
    if (i == refcount): 
    # When last ref molecule is reached
        sortoblarr = np.copy(oblarr)
        # Sort oblarr array
        sortoblarr.sort()
       
print ("\nQuery: ", Chem.MolToSmiles(query))
# Print query smile        
f = open("%s.dat" % output_name, "a") 
# Open output file and write query smile and fragments generated
f.write("\n\nQuery: {}\n{} fragments generated\n".format(Chem.MolToSmiles(query),num_frag))
f.close()

matcharr = []  
# list to store match molecule name
matchsmiarr = []
# list to store match molecule smile 
fragnumarr = np.zeros(4)
# array to store match molecule fragment number
refcount = 0
# Initialise ref count
for mol in range(len(refmollist)):
    fragnum = int(oblarr[refcount])
    name = reffraglist[refcount].name
    smile = reffraglist[refcount].smile
    # Store values of current molecule's name, smile, and fragment number
    
    for j in range(4): 
    # For the 4 largest fragment match
        idx = len(refmollist) - j - 1
        if (fragnum == sortoblarr[idx]): 
        # If refmol fragment number matches largest match fragment number
            matcharr.insert(j, name) 
            matchsmiarr.insert(j, smile)
            fragnumarr[j] = fragnum
            
    
    if (refcount == len(sortedsize) - 1):
        for i in range (4):
            print ("\nFrag match {}: {}\n {} with {} match fragments".format(i+1,matcharr[i], matchsmiarr[i], int(fragnumarr[i])))
            # Print 4 largest fragment matches
            f = open("%s.dat" % output_name, "a") 
            # Open output file and write 4 largest fragment matches
            f.write("\nFrag match {}: {}\n {} with {} match fragments\n".format(i+1,matcharr[i], matchsmiarr[i], int(fragnumarr[i])))
            f.close()
     
    refcount = refcount + 1
                         
                         ##Count number of C, N, and O##

queryatom_count = np.zeros(3)
diff =  np.zeros((len(refmollist),3))
# Count C, N, and O in query molecule
queryatom_count[0] = xyz.ele_count(query, "C")
queryatom_count[1] = xyz.ele_count(query, "N")
queryatom_count[2] = xyz.ele_count(query, "O")
totdiff = np.zeros(len(refmollist)) 
# Array for storing total C, N, and O number difference
sorttotdiff = np.zeros(len(refmollist))
# array for sorting totdiff

refatom_count = np.zeros((len(refmollist),3))
i = 0
for mol in ref:
    refatom_count[i][0] = xyz.ele_count(mol, "C")
    refatom_count[i][1] = xyz.ele_count(mol, "N")
    refatom_count[i][2] = xyz.ele_count(mol, "O")
    # Count C, N, and O in current ref molecule
    for j in range(3):
        if (refatom_count[i][j] > queryatom_count[j]):
            diff[i][j] =  refatom_count[i][j] - queryatom_count[j]
        else:
            diff[i][j] =  queryatom_count[j] - refatom_count[i][j]
            # Calculate ref and query atom number difference
        totdiff[i] = int(totdiff[i] + diff[i][j])
        # Add up difference for C, N, and O
    i = i+1
    if (i == len(refmollist)): 
        sorttotdiff = np.copy(totdiff)
        sorttotdiff.sort()
        # When last molecule is reached, sort totdiff

i = 0
matcharr = [] 
diffarr = np.zeros(4) 
# list to store match molecule name
for mol in ref:
    name = mol.GetProp("_Name")    
    diff = totdiff[i]
    
    for j in range(4): 
    # For the 4 least atom number difference
        if (diff == sorttotdiff[j]): 
        # If atom diff matches smallest atom diff
            matcharr.insert(j, name) 
            diffarr[j] = diff
        
    if (i == len(refmollist) - 1):
        for i in range (4):
            print ("\nAtom match {}: {}\n with {} difference".format(i+1,matcharr[i], int(diffarr[i])))
            # Print 4 smallest atom diff
            f = open("%s.dat" % output_name, "a") 
            # Open output file and write 4 smallest atom diff
            f.write("\nAtom match {}: {}\n with {} difference\n".format(i+1,matcharr[i], int(diffarr[i])))
            f.close()
     
    i = i + 1
        
        



    
     