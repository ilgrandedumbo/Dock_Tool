# ==============================================================================
# TL_Functions_modify_2.py - Part of Dock Tool (dock_tool)
# Main License: Creative Commons Attribution-NonCommercial-ShareAlike 4.0 (CC BY-NC-SA 4.0)
# Copyright (c) 2026 Papini, Sardelli.
# ------------------------------------------------------------------------------
# THIRD-PARTY CODE LICENSE NOTICE:
# Portions adapted from DOCK3.7 (Copyright 2025 The Regents of the University of California)
# distributed under the BSD-3-Clause License. See full disclaimer in this header.
# ==============================================================================
# THIRD-PARTY CODE EMBEDDED LICENSE NOTICE:
#
# Portions of this file (specifically molecule/DB2 parsing logic and suppliers)
# are adapted from code originally distributed under the following license:
#
# Copyright 2025 The Regents of the University of California
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, 
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, 
#    this list of conditions and the following disclaimer in the documentation 
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors 
#    may be used to endorse or promote products derived from this software 
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE.
# ==============================================================================

import xml.etree.ElementTree as ET #For reading XML files
from rdkit import Chem #The module with the RDKit functions we will need
from rdkit.Chem import AllChem #For 3D coordinate generation
import os #For the function to read in the molecules from the .mol2 file

# Import the XML file, using this as a guide:
# https://docs.python.org/3/library/xml.etree.elementtree.html
tree = ET.parse("TL_2.1_VERSION_6.xml")
root = tree.getroot()

# We will turn the procedure for estimating torsion strain energy in
# TL_Lookup_Test.py into a function that we can call on every molecule
# in the list, which will return an object of the following class:
class TP_list(object):
   def __init__(self, indices, exceptions, smarts, hc, methods, optimal_angles, hist_energy, actual_energy, mu):
       self.indices = indices #List of lists of indices
       self.exceptions = exceptions  # List of exception counts per torsion pattern
       self.smarts = smarts #List of SMARTS strings
       self.hc = hc #List of strings for hierarchy class
       self.methods = methods #List of strings of energy-estimation methods
       self.optimal_angles = optimal_angles #List of optimal angles
       self.TP_indices = [j for j in range(len(indices))]
       self.hist_energy = hist_energy
       self.actual_energy = actual_energy
       self.mu = mu
       # Indices of the torsion patterns, not the atoms!


   # Getters and setters:
   def get_indices(self):
       return(self.indices)
   def set_indices(self, inds):
       self.indices = inds
   def get_smarts(self):
       return(self.smarts)
   def set_smarts(self, sms):
       self.smarts = sms
   def get_hc(self):
       return(self.hc)
   def set_hc(self, hcs):
       self.hc = hcs
   def get_methods(self):
       return(self.methods)
   def set_methods(self, meths):
       self.methods = meths
   def get_optimal_angles(self):
       return(self.optimal_angles)
   def set_optimal_angles(self, angs):
    self.optimal_angles = angs
   def get_TP_indices(self): #Don't need a setter for this
       return(self.TP_indices)
   def get_hist_energy(self):
       return(self.hist_energy)
   def set_hist_energy(self, energy):
       self.hist_energy = energy
   def get_actual_energy(self):
       return(self.actual_energy)
   def set_actual_energy(self, energy):
       self.actual_energy = energy
   def get_max_energy(self):
       return(self.mu)
   def set_mu(self, energy):
       self.mu = energy


   # A method to return the information for a subset of the torsion patterns
   def get_TPs(self, inds = None):
       # The parameter inds is a list of indices for the torsion patterns,
       # not the atoms! We default to using all of the torsion patterns
       if inds == None:
           inds = [j for j in range(len(self.indices))]
       # Create a list of torsion pattern info to be returned:
       tps = [] #Initialize
       for j in inds:
           tps.append([self.TP_indices[j], self.exceptions[j], self.indices[j],
           self.smarts[j], self.hc[j], self.methods[j], self.optimal_angles[j], self.hist_energy[j], self.actual_energy[j], self.mu[j]])
       return(tps)


def smiMolSupplier(file_path):
    names = []
    ms = {}
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    used_names = set()

    with open(file_path, "r", encoding="utf-8") as f:
        for idx, raw_line in enumerate(f, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()
            smiles = parts[0]
            name = parts[1] if len(parts) > 1 else base_name

            if name in used_names:
                suffix = 2
                candidate_name = f"{name}_{suffix}"
                while candidate_name in used_names:
                    suffix += 1
                    candidate_name = f"{name}_{suffix}"
                name = candidate_name

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"[!] Warning: invalid SMILES at line {idx}: {smiles}")
                continue

            mol = Chem.AddHs(mol)

            params = AllChem.ETKDGv3()
            params.randomSeed = 0xF00D

            if AllChem.EmbedMolecule(mol, params) != 0:
                print(f"[!] Warning: 3D embedding failed for {name}")
                continue

            try:
                AllChem.UFFOptimizeMolecule(mol)
            except Exception:
                pass

            mol.SetProp("_Name", name)
            names.append(name)
            ms[name] = mol
            used_names.add(name)

    return names, ms


def sdfMolSupplier(file_path):
    names = []
    ms = {}

    supplier = Chem.SDMolSupplier(file_path, removeHs=False, sanitize=False)
    for idx, mol in enumerate(supplier, start=1):
        if mol is None:
            print(f"[!] Warning: failed to read SDF molecule at index {idx} in {file_path}")
            continue

        if mol.HasProp("_Name") and mol.GetProp("_Name").strip():
            name = mol.GetProp("_Name").strip()
        else:
            name = f"sdfNumber{idx}"

        mol.SetProp("_Name", name)
        names.append(name)
        ms[name] = mol

    return names, ms


def Mol2MolSupplier (file = None):
   names = [] #Make a list to hold the molecule names
   mols = {} #Make a dictionary
   with open(file, 'r') as f:
       fileend = os.fstat(f.fileno()).st_size
       count = 0
       line = f.readline()
       while not f.tell() == fileend:
           if line.startswith("#") or line == '\n':
               line = f.readline()
           if line.startswith("@<TRIPOS>MOLECULE"):
               count += 1
               mol = []
               mol.append(line)
               line = f.readline()
               if line != "\n" and line.split()[0].strip() not in names:
                   name = line.split()[0].strip()
               else:
                   name = "mol2Number" + str(count)


               while not line.startswith("@<TRIPOS>MOLECULE"):
                   mol.append(line)
                   line = f.readline()


                   if f.tell() == fileend:
                       mol.append(line)
                       break
               block = ",".join(mol).replace(',','')
               m = Chem.rdmolfiles.MolFromMol2Block(block, sanitize=False, removeHs = False)
               if m is not None:
                   m.SetProp("_Mol2Block", block)
                   m.SetProp("_Mol2Name", name)
               names.append(name)
               mols[name] = m
   return(names, mols)


# Here is an updated version to use with the output "file" string buffer object
# created by the db2_file_like function in the db2_to_mol2.py script
def db2MolSupplier(file):
   names = [] #Make a list to hold the molecule names
   mols = {} #Make a dictionary
   with file as f: #file is already opened as a string buffer
       bufferend = len(f.getvalue())
       count = 0
       line = f.readline()
       while not f.tell() == bufferend:
           if line.startswith("#") or line == '\n':
               line = f.readline()
           if line.startswith("@<TRIPOS>MOLECULE"):
               count += 1
               mol = []
               mol.append(line)
               line = f.readline()


               name = "db2Number" + str(count)


               while not line.startswith("@<TRIPOS>MOLECULE"):
                   mol.append(line)
                   line = f.readline()


                   if f.tell() == bufferend:
                       mol.append(line)
                       break
               block = ",".join(mol).replace(',','')
               m = Chem.rdmolfiles.MolFromMol2Block(block, sanitize=False, removeHs = False)
               if m is not None:
                   m.SetProp("_Mol2Block", block)
                   m.SetProp("_Mol2Name", name)
               names.append(name)
               mols[name] = m
   return(names, mols)



# This function will allow us to do the matching for each torsion rule
def tp_match(tp, hc, j, mol, bi):
   # tp is a torsion pattern, hc is the type of hierarchyClass ("general" or "specific"),
   # j is the current value for i, mol is the RDKit molecule,
   # bi is the bond_info list, threshold is the max energy allowed.

    exception = 0

    smarts = tp.get("smarts")
  
   # Create the histograms for energy estimates and bounds of confidence
   # intervals, if available
    pattern = Chem.MolFromSmarts(smarts)

    hist_E = []

    if tp.get("method") == "exact":
        for bin in tp.find("histogram_converted").findall("bin"):
            hist_E.append(float(bin.get("energy")))
        threshold = sum(hist_E) / len(hist_E)
    matches = mol.GetSubstructMatches(pattern)
   # A list of lists
    for match in matches: #For each match
       # Some of the SMARTS for the torion patterns actually have 5 atoms.
       # We need to ingore these
       if len(match) > 4:
           exception +=1
           continue #Go to the next match
       if mol.GetAtomWithIdx(match[0]).GetSymbol()=='H' or mol.GetAtomWithIdx(match[3]).GetSymbol()=='H':
           continue
       # Changed next line from "TP.get" to "tp.get"
       if tp.get("method") == "exact": #If using the exact method
            optimal_angles = [] # Initialize for the discretized angles
            actual_energy = [] # Initialize for the actual energy of this match
# Check if this minimum is within our acceptable energy threshold and min method
            N = len(hist_E)
            for i in range(N):
                prev = (i-1)%N
                next = (i+1)%N
                if hist_E[i] == 0 or (hist_E[i]<=hist_E[prev] and hist_E[i]<=hist_E[next]):
                    if hist_E[i] <= threshold:
                        center = -170.0 + i * 10.0
                        if [center-6, center+6] not in optimal_angles:
                            optimal_angles.append([center - 6.0, center + 6.0])
            bi.append(
           [
           list(match), #Convert tuple to list
           exception,
           smarts,
           hc, #"general" or "specific"
           "exact",
           optimal_angles,
           j, 
           False, #We will take this out when we create the final object
           hist_E,
           actual_energy,
           threshold
           ]
           )
       else: #If using the approximate method
           optimal_angles = [] # Initialize for the discretized angles
           actual_energy = [1.0, 0]
           threshold = None
           hist_E = None
           for angle in tp.find("angleList").findall("angle"):
               theta_0 = float(angle.get("theta_0")) #Peak location
               tolerance = float(angle.get("tolerance1")) #Tolerance around the peak
               optimal_angles.append([theta_0 - tolerance, theta_0 + tolerance])


           bi.append(
           [
           list(match), #Convert tuple to list
           exception,
           smarts,
           hc, #"general" or "specific"
           "approximate",
           optimal_angles,
           j,
           False, #We will take this out when we create the final object
           hist_E,
           actual_energy, 
           threshold 
           ]
           )




# This most general function will automate what we did in TL_Lookup_Test.py
def TL_lookup(mol): #mol is read in from the .mol2 file


   # List of lists of atom coordinates. Luckily RDKit starts indexing at 0
   bond_info = []
   # Initialize an empty list that will hold the information for each bond
   i = 0 #Initialize count of torsion rules
   # Loop over all of the specific hierarchy classes
   for HC in root.findall("hierarchyClass"):
       if HC.get("name") != "GG": #Not the general class
           for TP in HC.iter("torsionRule"): #Loop over each torsion rules
               tp_match(TP, "specific", i, mol, bond_info)
               i += 1 #Increase the count for the torsion rule


   # Now for the general method:
   for TP in root.find("hierarchyClass[@name='GG']").iter("torsionRule"):
       tp_match(TP, "general", i, mol, bond_info)
       i += 1 #Increase the count for the torsion rule

   if len(bond_info) == 0:
       return TP_list([], [], [], [], [], [], [], [], [])
       
   # Now that we have all of the torsion patterns, we need to be able to find
   # duplicates. The first such way is if the entire pattern is reversed. We
   # can fix this problem by making sure that all of the lists of indeces have
   # the second index (the first in the bond of interest) lower than the third
   # index (the second in the bond of interest)
   for  bond in bond_info: #Loop over every bond
       if bond[0][1] > bond[0][2]:
           bond[0].reverse() #Reverse this list
           bond[7] = True #Mark that we reversed this bond's indeces
       else:
           bond[7] = False #Mark that we did not reverse this bond's
           # indeces. We will remove this marking later


   # Next we condense the bond_info by the lists of 4 atoms defining the bonds.
   # We will pick the entry of bond_info that has the lowest value for i for
   # each match, since the torsion rules in the Torsion Library are arranged
   # (within each hierarchy class or hierarchy subclass) in decreasing
   # specificity, and we loop over all of the specific hierarchy classes
   # before the general one
   bond_info_red = [bond_info[0]] #Initialize a list for the reduced bond info
   # This reduced list needs at least one element for checking subelements
   for j in range(1, len(bond_info)):
       # Skip the first bond, which is already in the reduced list
       atom_0 = bond_info[j][0][0] #First atom index
       atom_1 = bond_info[j][0][1] #Second atom index
       atom_2 = bond_info[j][0][2] #Third bond index
       atom_3 = bond_info[j][0][3] #Fourth bond index


       unmatched = True #Initialize not finding a match
       for k in range(len(bond_info_red)):
           # Check against everything in the growing reduced list
           if bond_info_red[k][0][0] == atom_0 \
           and bond_info_red[k][0][1] == atom_1 \
           and bond_info_red[k][0][2] == atom_2 \
           and bond_info_red[k][0][3] == atom_3:
               # If there is a match in ALL of the atom indeces
               unmatched = False
               if bond_info[j][6] < bond_info_red[k][6]:
                   # Index 6 gives the torsion rule number j.
                   # If the new bond has a lower value, then we use it
                   # to replace the current one
                   bond_info_red[k] = bond_info[j]
                   break
                   # No need to continue looking for matches, since there
                   # should be no more than 1


       if unmatched: #If no match
           bond_info_red.append(bond_info[j]) #Append the current bond


   central_bonds = {}
   for j in range(len(bond_info_red)):
       bond = bond_info_red[j]
       atom_1 = bond[0][1]
       atom_2 = bond[0][2]
       key = (atom_1, atom_2) 

       if key not in central_bonds:
           central_bonds[key] = {
               'best_bond': bond,
               'all_matches': [bond[0]] 
           }
       else:
           current_best = central_bonds[key]['best_bond']
           
           is_better = False
           if bond[3][0] > current_best[3][0]: # 's' > 'g' (specific > general)
               is_better = True
           elif bond[3][0] == current_best[3][0] and bond[6] < current_best[6]:
               is_better = True
               
           is_same_rule = (bond[3] == current_best[3] and bond[6] == current_best[6])

           if is_better:
               central_bonds[key] = {
                   'best_bond': bond,
                   'all_matches': [bond[0]]
               }
           elif is_same_rule:
               central_bonds[key]['all_matches'].append(bond[0])

   b_i_r = []
   for key, data in central_bonds.items():
       best_bond = data['best_bond']
       matches = data['all_matches']
       
       if best_bond[7]: 
           for m in matches:
               m.reverse()
               
       best_bond[0] = matches 
       b_i_r.append(best_bond)

   return(
       TP_list(
           [bond[0] for bond in b_i_r], #List of lists of indeces for each bond
           [bond[1] for bond in b_i_r], #List of exception counts for each bond
           [bond[2] for bond in b_i_r], #List of SMARTS for each bond
           [bond[3] for bond in b_i_r], #List of hc class 
           [bond[4] for bond in b_i_r], #List of method for each bond 
           [bond[5] for bond in b_i_r], #List of optimal angles for each bond
           [bond[8] for bond in b_i_r], #hist_E
           [bond[9] for bond in b_i_r], #actual_energy
           [bond[10] for bond in b_i_r]  #mu
       ) 
    )
