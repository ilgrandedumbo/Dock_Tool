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
import csv
import os
import re
import torsion_pattern_analyzer as TPA
from statistics import pstdev
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem

def sanitize_filename(name):
    """Replace problematic characters for filesystem with underscores."""
    sanitized = re.sub(r'[=()[\]@#:/\\?\*<>|"]', '_', name)
    sanitized = re.sub(r'_+', '_', sanitized)
    sanitized = sanitized.strip('_')
    return sanitized

def calculate_sigma_mu(hist_energy, mu):
    if hist_energy is None:
        return "NA"
    return pstdev(hist_energy) / mu

def energy_sort_key(value):
    if value == "NA":
        return float("inf")
    return value

def generate_csv_reports(input_file, results_list, molecule_dict=None, create_matrix=True, disable_progress=False):
    # Extract directory and base name from input_file
    output_dir = os.path.dirname(input_file)
    base_path = os.path.basename(input_file)
    
    # Determine output name based on extension
    if base_path.endswith(".mol2"):
        output_name = base_path[:-5]
    elif base_path.endswith(".db2"):
        output_name = base_path[:-4]
    elif base_path.endswith(".db2.gz"):
        output_name = base_path[:-7]
    elif base_path.endswith(".smi"):
        output_name = base_path[:-4]
    elif base_path.endswith(".sdf"):
        output_name = base_path[:-4]
    else:
        output_name = os.path.splitext(base_path)[0]

    # Create output directory if it doesn't exist
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Sanitize only the filename (not the path) for the filesystem
    safe_output_name = sanitize_filename(output_name)
    raw_file_name = os.path.join(output_dir, safe_output_name + "_Torsion_Strain_raw.csv")
    clean_file_name = os.path.join(output_dir, safe_output_name + "_Torsion_Strain_clean.csv")
    matrix_file_name = os.path.join(output_dir, safe_output_name + "_Torsion_Strain_matrix.csv")
    energy_file_name = os.path.join(output_dir, safe_output_name + "_Torsion_Strain_energy.csv")

    optimal_angles_count = {} 
    no_optimal_angle_count = {}
    na_count = 0

    prepared_results = []

    for raw_result in results_list:
        if isinstance(raw_result, dict):
            name = raw_result.get("name", "unknown")
            analysis_row = raw_result.get("analysis", [])
            bond_data_list = analysis_row[1:] if analysis_row else []
            molecule_entry = raw_result.get("molecule")
        else:
            name = raw_result[0]
            bond_data_list = raw_result[1:]
            molecule_entry = molecule_dict.get(name) if molecule_dict else None

        try:
            analysis_row = analysis_row if isinstance(raw_result, dict) else raw_result

            if len(analysis_row) > 1 and analysis_row[1] == "NA":
                prepared_results.append({
                    "name": name,
                    "analysis_row": analysis_row,
                    "bond_data_list": bond_data_list,
                    "molecule_entry": molecule_entry,
                    "global_energy": "NA",
                })
                continue

            global_energy = "NA"
            vdw, elec, torsion, bond_term, angle_term = "NA", "NA", "NA", "NA", "NA"

            if molecule_entry is None and molecule_dict and name in molecule_dict:
                molecule_entry = molecule_dict[name]

            if molecule_entry is not None:
                try:
                    rdkit_mol = molecule_entry
                    mol_hs = Chem.AddHs(rdkit_mol)
                    mp = AllChem.MMFFGetMoleculeProperties(mol_hs)
                    ff = AllChem.MMFFGetMoleculeForceField(mol_hs, mp)

                    if ff is not None:
                        global_energy = ff.CalcEnergy()

                        mp.SetMMFFVdWTerm(False)
                        vdw = global_energy - AllChem.MMFFGetMoleculeForceField(mol_hs, mp).CalcEnergy()
                        mp.SetMMFFVdWTerm(True)

                        mp.SetMMFFEleTerm(False)
                        elec = global_energy - AllChem.MMFFGetMoleculeForceField(mol_hs, mp).CalcEnergy()
                        mp.SetMMFFEleTerm(True)

                        mp.SetMMFFTorsionTerm(False)
                        torsion = global_energy - AllChem.MMFFGetMoleculeForceField(mol_hs, mp).CalcEnergy()
                        mp.SetMMFFTorsionTerm(True)

                        mp.SetMMFFBondTerm(False)
                        bond_term = global_energy - AllChem.MMFFGetMoleculeForceField(mol_hs, mp).CalcEnergy()
                        mp.SetMMFFBondTerm(True)

                        mp.SetMMFFAngleTerm(False)
                        angle_term = global_energy - AllChem.MMFFGetMoleculeForceField(mol_hs, mp).CalcEnergy()
                        mp.SetMMFFAngleTerm(True)
                except Exception as ef:
                    print(f"Warning: Impossibile calcolare FF per {name}: {ef}")

            prepared_results.append({
                "name": name,
                "analysis_row": analysis_row,
                "bond_data_list": bond_data_list,
                "molecule_entry": molecule_entry,
                "global_energy": global_energy,
                "vdw": vdw,
                "elec": elec,
                "torsion": torsion,
                "bond_term": bond_term,
                "angle_term": angle_term,
            })
        except Exception as e:
            na_count += 1
            print(f"Error processing {name}: {e}")
            prepared_results.append({
                "name": name,
                "analysis_row": [name, "NA"],
                "bond_data_list": [],
                "molecule_entry": None,
                "global_energy": "NA",
                "vdw": "NA",
                "elec": "NA",
                "torsion": "NA",
                "bond_term": "NA",
                "angle_term": "NA",
            })

    prepared_results.sort(key=lambda item: energy_sort_key(item.get("global_energy", "NA")))

    # Open files for Raw and Clean data
    with open(raw_file_name, mode="w", newline='') as raw_file, \
         open(clean_file_name, mode="w", newline='') as clean_file, \
         open(energy_file_name, mode="w", newline='') as energy_file:
        
        raw_writer = csv.writer(raw_file)
        clean_writer = csv.writer(clean_file)
        energy_writer = csv.writer(energy_file)
        # Headers
        clean_writer.writerow(["Molecule Name", "Total Rotable Axes", "Rotable in Minimum", "Skipped 4+ Atom Matches"])
        raw_writer.writerow(["Molecule Name", "Central Axes", "Representative Smart", "Branch Atoms", "Branch Angle", "Optimal angle(s)", "Actual Energy", "Mu", "Sigma/Mu", "Is Optimal"])
        energy_writer.writerow(["Molecule Name", "GLOBAL_MMFF94", "Van_der_Waals", "Electrostatic", "Torsion_Strain", "Bond_Stretch", "Angle_Bend"])       

        for item in tqdm(prepared_results, desc="Processing analysis results", disable=disable_progress):
            name = item["name"]
            analysis_row = item["analysis_row"]
            bond_data_list = item["bond_data_list"]
            molecule_entry = item["molecule_entry"]
            global_energy = item.get("global_energy", "NA")
            vdw = item.get("vdw", "NA")
            elec = item.get("elec", "NA")
            torsion = item.get("torsion", "NA")
            bond_term = item.get("bond_term", "NA")
            angle_term = item.get("angle_term", "NA")

            if len(analysis_row) > 1 and analysis_row[1] == "NA":
                na_count += 1
                clean_writer.writerow([name, "NA", "NA", "NA"])
                raw_writer.writerow([name, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"])
                energy_writer.writerow([name, "NA", "NA", "NA", "NA", "NA", "NA"])
                continue

            energy_writer.writerow([name, global_energy, vdw, elec, torsion, bond_term, angle_term])

            try:
                for bond in bond_data_list:
                    smart = bond[3]
                    atoms = bond[2][0]
                    central_axis = bond[2][0][1:3]
                    actual_angle = bond[7][0]
                    optimal_angles = bond[6]
                    is_in_optimal = bond[8]
                    actual_energy = bond[9]
                    mu = bond[10]
                    hist_E = bond[11]
                    sigma_mu = calculate_sigma_mu(hist_E, mu)
                    
                    if create_matrix:
                        if smart not in optimal_angles_count:
                            optimal_angles_count[smart] = {}
                        if smart not in no_optimal_angle_count:
                            no_optimal_angle_count[smart] = {}
                            
                        if is_in_optimal:
                            matched_min = None
                            for ang_min in optimal_angles:
                                if TPA.angle_within_tol(actual_angle, ang_min):
                                    matched_min = tuple(ang_min)
                                    break
                            if matched_min is not None:
                                d = optimal_angles_count[smart]
                                if matched_min not in d:
                                    d[matched_min] = []
                                d[matched_min].append(sigma_mu)
                        else:
                            ang_key = float(int(actual_angle))
                            d = no_optimal_angle_count[smart]
                            if ang_key not in d:
                                d[ang_key] = []
                            d[ang_key].append(sigma_mu)
                            
                    centered_optimal_angles = []
                    if isinstance(optimal_angles, list):
                        for r in optimal_angles:
                            if isinstance(r, (list, tuple)) and len(r) == 2:
                                center = (r[0] + r[1]) / 2.0
                                centered_optimal_angles.append(round(center, 1))
                            else:
                                centered_optimal_angles.append(r)
                    else:
                        centered_optimal_angles = optimal_angles

                    raw_writer.writerow([name, central_axis, smart, atoms, actual_angle, centered_optimal_angles, actual_energy, mu, sigma_mu, is_in_optimal])
                
                tot_matches = len(bond_data_list)
                tot_at_min = sum(1 for b in bond_data_list if b[8] == True)
                tot_eccezioni = sum(b[1] for b in bond_data_list)

                clean_writer.writerow([name, tot_matches, tot_at_min, tot_eccezioni])

            except Exception as e:
                na_count += 1
                print(f"Error processing {name}: {e}")
                raw_writer.writerow([name, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"])
                clean_writer.writerow([name, "NA", "NA", "NA"])
                energy_writer.writerow([name, "NA", "NA", "NA", "NA", "NA", "NA"])

    # Create Matrix CSV only if create_matrix is True
    if create_matrix:
        with open(matrix_file_name, mode="w", newline='') as matrix_file:
            matrix_writer = csv.writer(matrix_file)
            matrix_writer.writerow(["Smart", "type", "angle(s) + frequency + sigma/mu"])
            
            for smart in optimal_angles_count.keys():
                try:
                    optimal_dict = optimal_angles_count.get(smart, {})
                    no_optimal_dict = no_optimal_angle_count.get(smart, {})

                    def process_data(data_dict, is_optimal=False):
                        sorted_items = sorted(data_dict.items(), key=lambda x: len(x[1]), reverse=True)
                        rows = []
                        for ang, sigma_mu_list in sorted_items:
                            count = len(sigma_mu_list)
                            total = sum(pair for pair in sigma_mu_list)
                            avg = total / count if count > 0 else float('inf')

                            
                            if is_optimal:
                                label = (ang[0] + ang[1]) / 2
                            else:
                                label = ang
                            
                            rows.append(f"{label:.1f}, {count}, {avg:.3f}")
                        return "; ".join(rows)

                    matrix_writer.writerow([smart, "optimal", process_data(optimal_dict, True)])
                    matrix_writer.writerow([smart, "no_optimal", process_data(no_optimal_dict, False)])
                except Exception as e:
                    print(f"Error in matrix for {smart}: {e}")
                    matrix_writer.writerow([smart, "optimal", "NA"])
                    matrix_writer.writerow([smart, "no_optimal", "NA"])

    # print(f"[*] Raw output -> {raw_file_name}")
    # print(f"[*] Clean output -> {clean_file_name}")
    # if create_matrix:
    #     print(f"[*] Matrix output -> {matrix_file_name}")
    
    # success_count = len(results_list) - na_count
    # print(f"[*] Analysis finished: {success_count} successful / {na_count} NA")