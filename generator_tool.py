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
import math
import itertools
import os
import sys
import random
import re
from tqdm import tqdm

import TL_Functions_modify_2 as TL
import torsion_pattern_analyzer as TPA

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    import rdkit.Chem as Chem
    from rdkit.Chem import AllChem

max_cap = 10000

def sanitize_filename(name):
    """Replace problematic characters for filesystem with underscores."""
    sanitized = re.sub(r'[=()[\]@#:/\\?\*<>|"]', '_', name)
    sanitized = re.sub(r'_+', '_', sanitized)
    sanitized = sanitized.strip('_')
    return sanitized

def calc_combination(bond_info):
    total_combination = 1
    for bond in bond_info:
        num_atom = len(bond[2])
        num_ang = len(bond[6])
        row_combination = num_ang * num_atom
        if row_combination > 0:
            total_combination *= row_combination
    return total_combination

def safe_mol_to_molblock(mol):
    """Convert an RDKit molecule to a MolBlock while handling aromaticity issues."""
    try:
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        mol.UpdatePropertyCache(strict=False)
    
    try:
        return Chem.MolToMolBlock(mol, kekulize=False)
    except Exception:
        try:
            return Chem.MolToMolBlock(mol)
        except Exception as e:
            raise RuntimeError("Unable to convert RDKit molecule to MolBlock due to aromaticity/kekulization errors") from e


def load_input_molecules(file_path):
    if file_path.endswith(".mol2"):
        return TL.Mol2MolSupplier(file_path)
    if file_path.endswith(".db2") or file_path.endswith(".db2.gz"):
        return TL.db2MolSupplier(file_path)
    if file_path.endswith(".sdf"):
        return TL.sdfMolSupplier(file_path)
    raise ValueError(f"Unsupported input format: {file_path}")


def extract_mol2_template(file_path, target_name):
    current_block = []
    current_name = None

    with open(file_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("@<TRIPOS>MOLECULE"):
                if current_block and current_name == target_name:
                    return current_block
                current_block = [line]
                current_name = None
                continue

            if current_block:
                current_block.append(line)
                if current_name is None and line.strip() and not line.startswith("@<TRIPOS>"):
                    current_name = line.strip().split()[0]

    if current_block and current_name == target_name:
        return current_block

    return None


def rewrite_mol2_coordinates(template_lines, mol):
    conformer = mol.GetConformer()
    atom_index = 0
    inside_atom_section = False
    output_lines = []

    for line in template_lines:
        stripped = line.strip()

        if line.startswith("@<TRIPOS>ATOM"):
            inside_atom_section = True
            output_lines.append(line)
            continue

        if inside_atom_section and line.startswith("@<TRIPOS>"):
            inside_atom_section = False

        if inside_atom_section and stripped:
            line_ending = "\n" if line.endswith("\n") else ""
            parts = line.split()
            if len(parts) >= 6:
                if atom_index >= conformer.GetNumAtoms():
                    raise ValueError("Mol2 template has more atom lines than the generated conformer")

                position = conformer.GetAtomPosition(atom_index)
                parts[2] = f"{position.x: .4f}"
                parts[3] = f"{position.y: .4f}"
                parts[4] = f"{position.z: .4f}"
                line = " ".join(parts) + line_ending
                atom_index += 1

        output_lines.append(line)

    if atom_index != conformer.GetNumAtoms():
        raise ValueError("Mol2 template has fewer atom lines than the generated conformer")

    return output_lines


def rename_mol2_title(template_lines, generated_name):
    renamed_lines = []
    inside_molecule_section = False
    title_replaced = False

    for line in template_lines:
        if line.startswith("@<TRIPOS>MOLECULE"):
            inside_molecule_section = True
            renamed_lines.append(line)
            continue

        if inside_molecule_section and not title_replaced and line.strip() and not line.startswith("@<TRIPOS>"):
            line_ending = "\n" if line.endswith("\n") else ""
            line = f"{generated_name}{line_ending}"
            title_replaced = True

        if inside_molecule_section and line.startswith("@<TRIPOS>"):
            inside_molecule_section = False

        renamed_lines.append(line)

    return renamed_lines


def write_generated_structure(mol, output_path, mol_name, source_file_path=None, generated_name=None):
    output_name = generated_name or mol_name

    if source_file_path and source_file_path.lower().endswith(".mol2"):
        template_text = mol.GetProp("_Mol2Block") if mol.HasProp("_Mol2Block") else None
        if template_text:
            template_lines = template_text.splitlines(keepends=True)
        else:
            template_lines = extract_mol2_template(source_file_path, mol_name)
        if template_lines is None:
            raise RuntimeError(f"Unable to locate Mol2 block for {mol_name} in {source_file_path}")
        renamed_lines = rename_mol2_title(template_lines, output_name)
        rewritten_lines = rewrite_mol2_coordinates(renamed_lines, mol)
        with open(output_path, "w", encoding="utf-8", newline="") as handle:
            handle.writelines(rewritten_lines)
        return

    sdf_mol = Chem.Mol(mol)
    sdf_mol.SetProp("_Name", output_name)
    writer = Chem.SDWriter(output_path)
    try:
        writer.write(sdf_mol)
    finally:
        writer.close()

def mol2_generator(rdkit_mol, analysis, number_of_gen, mol_name, gacha_mode, input_file_base=None, source_file_path=None, force=False, strain_mode=False):
    bond_data_list = analysis[1:]
    combination = calc_combination(bond_data_list)
    if number_of_gen <= 0 or number_of_gen > combination:
        number_of_gen = combination
        
    count = 0
    safe_mol_name = sanitize_filename(mol_name)
    if input_file_base:
        safe_input_base = sanitize_filename(os.path.splitext(input_file_base)[0])
        output_dir = os.path.join("generated_conformers", safe_input_base, safe_mol_name)
    else:
        output_dir = f"output_{safe_mol_name}"
    
    os.makedirs(output_dir, exist_ok=True)
    
    all_bond_options = []
    for bond in bond_data_list:
        atoms_list = bond[2]
        optimal_angles = bond[6]
        flattened_angles = [(r[0] + r[1]) / 2 for r in optimal_angles]
        bond_options = []
        for atoms in atoms_list:
            for angle in flattened_angles:
                bin_num = ((math.ceil(angle / 10.0) + 17) % 36)
                actual_energy = bond[11][bin_num] if bond[11] is not None else float('inf')
                bond_options.append((atoms, angle, actual_energy))
        all_bond_options.append(bond_options)

    max_candidate_cap = None if force else max_cap
    candidate_combos = []
    if gacha_mode:
        viewed_identifiers = set()
        pool_size = max(number_of_gen * 5, 2000)
        if max_candidate_cap is not None:
            pool_size = min(pool_size, max_candidate_cap)
        num_of_try = pool_size * 10
        i = 0
        while len(candidate_combos) < pool_size and i < num_of_try:
            i += 1
            random_conformer = [random.choice(options) for options in all_bond_options]
            identifier = tuple((tuple(atoms), angle) for atoms, angle, _ in random_conformer)
            if identifier not in viewed_identifiers:
                viewed_identifiers.add(identifier)
                candidate_combos.append(random_conformer)
    else:
        max_combinations = combination if max_candidate_cap is None else min(combination, max_candidate_cap)
        combo_generator = itertools.product(*all_bond_options)
        candidate_combos = list(itertools.islice(combo_generator, max_combinations))

    # Added H atoms to ensure proper MMFF94 energy calculations, and precompute molecule properties for efficiency
    mol_with_hs = Chem.AddHs(rdkit_mol)
    mp = AllChem.MMFFGetMoleculeProperties(mol_with_hs)
    
    # Create a list to hold scored candidates for sorting and selection
    scored_candidates = []
    desc_text = f"Screening Strain per {mol_name}" if strain_mode else f"Screening MMFF94 in RAM per {mol_name}"
    for combo in tqdm(candidate_combos, desc=desc_text):
        mol_temp = Chem.Mol(mol_with_hs)
        conf = mol_temp.GetConformer()
        
        hist_energy_sum = 0
        conformer_geo = []
        for atoms, angle, h_energy in combo:
            Chem.rdMolTransforms.SetDihedralDeg(conf, atoms[0], atoms[1], atoms[2], atoms[3], float(angle))
            hist_energy_sum += h_energy
            conformer_geo.append((atoms, angle))
            global_mmff_energy = None
        if not strain_mode:   
            ff = AllChem.MMFFGetMoleculeForceField(mol_temp, mp)
            if ff is None:
                continue
            global_mmff_energy = ff.CalcEnergy()
        scored_candidates.append({
            'conformer_geo': conformer_geo,
            'hist_energy_sum': hist_energy_sum,
            'global_energy': global_mmff_energy,
            'mol_obj': mol_temp
        })
        
    if not strain_mode:
        scored_candidates.sort(key=lambda x: x['global_energy'])
    else:
        scored_candidates.sort(key=lambda x: x['hist_energy_sum'])
    winners = scored_candidates[:number_of_gen]
    
    if not force and not gacha_mode and len(candidate_combos) >= max_cap and number_of_gen >= max_cap:
        print(f"[!] Reached internal limit of {max_cap} combinations for {mol_name}. Generation stopped early.")

    log_file_path = os.path.join(output_dir, f"{safe_mol_name}.txt")
    if os.path.exists(log_file_path):
        os.remove(log_file_path)
        
    for idx, item in enumerate(tqdm(winners, desc=f"Writing top conformers for {mol_name}")):
        conformer_geo = item['conformer_geo']
        global_energy = item['global_energy']
        hist_energy_sum = item['hist_energy_sum']
        mol_rdkit = item['mol_obj']
        
        if strain_mode:
            mp_local = AllChem.MMFFGetMoleculeProperties(mol_rdkit)
            if mp_local is not None:
                ff_local = AllChem.MMFFGetMoleculeForceField(mol_rdkit, mp_local)
                if ff_local is not None:
                    global_energy = ff_local.CalcEnergy()
                    item['global_energy'] = global_energy
        
        # Decompose global energy into components by selectively disabling terms in the MMFF force field, this for a problem with the MMFFGetMoleculeForceField method that doesn't provide direct access to energy components
        mp_winner = AllChem.MMFFGetMoleculeProperties(mol_rdkit)
        
        mp_winner.SetMMFFVdWTerm(False)
        vdw = global_energy - AllChem.MMFFGetMoleculeForceField(mol_rdkit, mp_winner).CalcEnergy()
        mp_winner.SetMMFFVdWTerm(True)
        
        mp_winner.SetMMFFEleTerm(False)
        elec = global_energy - AllChem.MMFFGetMoleculeForceField(mol_rdkit, mp_winner).CalcEnergy()
        mp_winner.SetMMFFEleTerm(True)
        
        mp_winner.SetMMFFTorsionTerm(False)
        torsion = global_energy - AllChem.MMFFGetMoleculeForceField(mol_rdkit, mp_winner).CalcEnergy()
        mp_winner.SetMMFFTorsionTerm(True)
        
        mp_winner.SetMMFFBondTerm(False)
        bond = global_energy - AllChem.MMFFGetMoleculeForceField(mol_rdkit, mp_winner).CalcEnergy()
        mp_winner.SetMMFFBondTerm(True)
        
        mp_winner.SetMMFFAngleTerm(False)
        angle = global_energy - AllChem.MMFFGetMoleculeForceField(mol_rdkit, mp_winner).CalcEnergy()
        mp_winner.SetMMFFAngleTerm(True)
        
        generation_tags = []
        if gacha_mode:
            generation_tags.append("[Gacha]")
        if strain_mode:
            generation_tags.append("[strain]")
        generation_tag = "".join(generation_tags) + f"_gen_{idx+1}"
        output_label = f"{mol_name}{generation_tag}"
        safe_output_label = sanitize_filename(output_label)
        geo_details = ", ".join([f"Atoms_{i+1}: {a} Ang_{i+1}: {ang:.1f}°" for i, (a, ang) in enumerate(conformer_geo)])
        energy_details = f"| Hist_Energy_Sum: {hist_energy_sum:.4f} | GLOBAL_MMFF94: {global_energy:.4f} [VDW: {vdw:.4f}, ELEC: {elec:.4f}, TORSION_FF: {torsion:.4f}, BOND_BEND: {bond+angle:.4f}]"
        
        with open(log_file_path, "a", encoding="utf-8") as f:
            f.write(f"{output_label}: {geo_details} {energy_details}\n")
            
        if source_file_path and source_file_path.lower().endswith(".mol2"):
            output_extension = ".mol2"
        else:
            output_extension = ".sdf"

        filename = os.path.join(output_dir, f"{safe_output_label}{output_extension}")
        write_generated_structure(
            mol_rdkit,
            filename,
            mol_name,
            source_file_path=source_file_path,
            generated_name=output_label,
        )
        count += 1
        
    return count

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit()

    input_file = sys.argv[1]
    names, ms = load_input_molecules(input_file)
        
    try:
        n_gen = int(input("\nHow much generation? press 'Enter' for all combinations "))
    except ValueError:
        n_gen = -2

    for name in names:
        mol = ms[name]
        if mol is None: continue
            
        print(f"\n--- Analisi molecola: {name} ---")
        analysis_result = TPA.Torsion_Strain(name, mol, {name: mol.GetConformer().GetPositions()})
            
        if analysis_result[1] == "NA":
            print(f"Skip {name}: Analysis.")
            continue
            
        num_created = mol2_generator(mol, analysis_result, n_gen, name, gacha_mode=False, input_file_base=os.path.basename(input_file), source_file_path=input_file)
        print(f"Completed: {num_created} conformation created for {name}")