# ==============================================================================
# dock_tool.py - Main Orchestrator
# Part of the Docking Simplifier Tool
#
# Copyright (c) 2026 Niccolò, Papini, Sardelli.
# Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
# International License (CC BY-NC-SA 4.0).
#
# Free for Academic and Research use. Commercial and SaaS usage is STRICTLY PROHIBITED.
# ==============================================================================
import os
from concurrent.futures import ProcessPoolExecutor
from sys import argv #For passing commandline arguments
import csv #For writing a .csv file
import TL_Functions_modify_2 as TL
#import db2_to_mol2 as converter
import numpy as np
import math
from math import atan2, pi, sqrt, ceil
# We will import the other required modules when we run TL_Functions.py

# number of thread, core used by multi-threading
num_workers = max(1, os.cpu_count() // 2)
def main():
    script, input = argv #Input from commandline argument
    if input[-5:] == ".mol2": #If fed in a .mol2 file
        names, ms = TL.Mol2MolSupplier(input)
        # A list of names and a dictionary of mol2 objects with the names as keys
        # Make sure every mol2 object in the file has a different name, or not all
        # of them will end up in the final dictionary!
        output_name = input[:len(input)-5] #Everything but the ".mol2"
    elif input[-4:] == ".db2" or input[-7:] == ".db2.gz": #If fed in a .db2 file
        file = converter.db2_file_like(input)
        # Creates a string buffer object called "file" that looks like a .mol2 file
        names, ms = TL.db2MolSupplier(file)
        output_name = input[:len(input)-4] #Everything but the ".db2"
    else:
        print("Error. Please pass a .mol2 or .db2 file.")
        exit()

    print(str(len(names)) + " molecules finished reading. Calculating strain energy...")

    positions_dict = {}
    for name, mol in ms.items():
        if mol is None or mol.GetNumConformers() == 0:
            positions_dict[name] = None
            continue
        positions = mol.GetConformer().GetPositions()
        positions_dict[name] = positions

    tasks = [(name, ms[name], positions_dict) for name in names]

    # multithread
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        elaborate_mol = list(executor.map(process_mol, tasks))
    print(f"{len(elaborate_mol)} molecules processed.")

def dist_to_ang(a, target):
    return abs((a - target + 180) % 360 - 180)

def calc_min_dist(angle, angles_min):
    min_dist = 180.0
    for low, upper in angles_min:
        if angle_within_tol(angle, [low, upper]):
            return 0.0
        d_low = dist_to_ang(angle, low)
        d_upper = dist_to_ang(angle, upper)
        min_dist = min(min_dist, d_low, d_upper)
    return min_dist

def angle_within_tol(angle, angles_min):
    # Normalize angle in [-180, 180]
    angle = float(angle)
    a = ((angle + 180) % 360) - 180
    low, upper = angles_min
    low = ((low + 180) % 360) - 180
    upper = ((upper + 180) % 360) - 180
    if low <= upper:
        if low <= a <= upper:
            return True
    else: 
        if a >= low or a <= upper:
            return True
    return False

def unit(a):
    return(a / sqrt(np.dot(a,a)))

def dihedral(a_1, a_2, a_3, a_4):
    b_1 = a_2 - a_1
    b_2 = a_3 - a_2
    b_3 = a_4 - a_3

    n_1 = unit(np.cross(b_1, b_2))
    n_2 = unit(np.cross(b_2, b_3))
    m = unit(np.cross(n_1, b_2))

    x = np.dot(n_1, n_2)
    y = np.dot(m, n_2)
    return -atan2(y, x) * 180 / pi

def Torsion_Strain(name, mol, positions_dict=None):
    if mol is not None: #Check to make sure the molecule exists
        try:
            M = TL.TL_lookup(mol) #Create a TP_list function
            mol_info = M.get_TPs() #The molecule's information
            bond_info = []
            for item in sorted(mol_info, key=lambda l: l[0], reverse=True):
                h = item[0]
                exceptions = item[1]
                atoms_matrix = item[2]
                smarts = item[3]
                tag1 = item[4]
                tag2 = item[5]
                ang_mins = item[6]
                hist_E = item[7]
                actual_energy = item[8]
                mu = item[9]
                
                if positions_dict[name] is not None:
                    pos = positions_dict[name]
                    pattern_results = []
                    for match_atoms in atoms_matrix:
                        i1, i2, i3, i4 = match_atoms
                        angle = dihedral(pos[i1], pos[i2], pos[i3], pos[i4])
                        below_threshold = any(angle_within_tol(angle, ang_min) for ang_min in ang_mins)
                        dist = calc_min_dist(angle, ang_mins)

                        pattern_results.append({
                                'atoms': match_atoms, 'angle': angle,
                                'is_optimal': below_threshold, 'dist': dist
                        })
                    pattern_results.sort(key=lambda x: x['dist'])
                    sorted_atoms = [p['atoms'] for p in pattern_results]
                    sorted_angles = [p['angle'] for p in pattern_results]
                    sorted_optimals = [p['is_optimal'] for p in pattern_results]
                    final_optimal = any(sorted_optimals)
                    theta = sorted_angles[0]
                    if hist_E is not None:
                        bin_num = ((math.ceil(theta / 10.0) + 17) % 36)
                        # Calcolo energia precisa al decimale
                        interpolated_energy = (hist_E[bin_num ] - hist_E[(bin_num + 35) % 36]) / 10.0 * (theta - (bin_num - 17) * 10) + hist_E[bin_num]
                        
                        closest_bin = int(round(((theta + 175) // 10)) % 36) 
                        actual_energy = (interpolated_energy, closest_bin)
                    else:
                        actual_energy = (1.0, -1)
                    bond_info.append([h, exceptions, sorted_atoms, smarts, tag1, tag2, ang_mins, sorted_angles, final_optimal, actual_energy, mu, hist_E])
                else:
                    bond_info.append(item)
            # Unlist the list of lists into just a list of elements
            return [name] + bond_info
        except Exception as e:
            print(f"Error processing {name}: {e}")
            return [name] + ["NA"]


def process_mol(args):
    name, mol, positions_dict = args
    return Torsion_Strain(name, mol, positions_dict)

# initialize the array with elaborate mol

elaborate_mol = []

if __name__ == "__main__":
    main()
    


# END
