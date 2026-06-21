# ==============================================================================
# dock_tool.py - Main Orchestrator
# Part of the Dock Tool
#
# Copyright (c) 2026 Niccolò, Papini, Sardelli.
# Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
# International License (CC BY-NC-SA 4.0).
#
# Free for Academic and Research use. Commercial and SaaS usage is STRICTLY PROHIBITED.
# ==============================================================================

import argparse
import gzip
from io import StringIO
import sys
import os
import re
import tempfile
import TL_Functions_modify_2 as TL
import torsion_pattern_analyzer as TPA
import Torsion_Strain_Writer as TSW
import generator_tool as generator
from tqdm import tqdm
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    import rdkit.Chem as Chem
    from rdkit.Chem import AllChem

def sanitize_filename(name):
    """Replace problematic characters for filesystem with underscores."""
    # Replace characters problematic for file/folder names
    sanitized = re.sub(r'[=()[\]@#:/\\?\*<>|"]', '_', name)
    # Remove consecutive underscores
    sanitized = re.sub(r'_+', '_', sanitized)
    # Remove leading/trailing underscores
    sanitized = sanitized.strip('_')
    return sanitized


def natural_sort_key(value):
    parts = re.split(r'(\d+)', value)
    return [int(part) if part.isdigit() else part.lower() for part in parts]


def write_temp_sdf_from_molecules(file_path, names, molecules):
    base_name = sanitize_filename(os.path.splitext(os.path.basename(file_path))[0])
    fd, temp_sdf_path = tempfile.mkstemp(prefix=f"{base_name}_", suffix=".sdf")
    os.close(fd)

    writer = Chem.SDWriter(temp_sdf_path)
    try:
        for name in names:
            mol = molecules.get(name)
            if mol is None:
                continue

            sdf_mol = Chem.Mol(mol)
            sdf_mol.SetProp("_Name", name)
            writer.write(sdf_mol)
    finally:
        writer.close()

    return temp_sdf_path


def load_normalized_molecules(file_path):
    if file_path.endswith(".mol2"):
        names, molecules = TL.Mol2MolSupplier(file_path)
        return file_path, names, molecules, None

    if file_path.endswith(".sdf"):
        names, molecules = TL.sdfMolSupplier(file_path)
        return file_path, names, molecules, None

    if file_path.endswith(".smi"):
        names, molecules = TL.smiMolSupplier(file_path)
    elif file_path.endswith(".db2.gz"):
        with gzip.open(file_path, "rt", encoding="utf-8", errors="replace") as handle:
            db2_buffer = StringIO(handle.read())
        names, molecules = TL.db2MolSupplier(db2_buffer)
    elif file_path.endswith(".db2"):
        with open(file_path, "r", encoding="utf-8", errors="replace") as handle:
            db2_buffer = StringIO(handle.read())
        names, molecules = TL.db2MolSupplier(db2_buffer)
    else:
        raise ValueError(f"Unsupported input format: {file_path}")

    temp_sdf_path = write_temp_sdf_from_molecules(file_path, names, molecules)
    names, molecules = TL.sdfMolSupplier(temp_sdf_path)
    return temp_sdf_path, names, molecules, temp_sdf_path


def collect_input_files(paths):
    """Collect supported molecule files from paths, grouping by top-level input item."""
    SUPPORTED_EXTENSIONS = (".mol2", ".db2", ".db2.gz", ".smi", ".sdf")
    SKIP_DIRS = {"analysis_results", "generated_conformers", "__pycache__"}
    
    grouped_files = {}
    for path_item in paths:
        if os.path.isdir(path_item):
            print(f"[*] Searching in directory: {path_item}")
            collected = []
            for root, dirs, files in os.walk(path_item):
                dirs[:] = [
                    d for d in dirs
                    if d not in SKIP_DIRS and not d.startswith("output_") and not d.startswith("csv_")
                ]
                for file_name in sorted(files, key=natural_sort_key):
                    if file_name.endswith(SUPPORTED_EXTENSIONS):
                        collected.append(os.path.join(root, file_name))
            grouped_files[path_item] = sorted(
                collected,
                key=lambda file_path: natural_sort_key(os.path.relpath(file_path, path_item))
            )
        elif os.path.isfile(path_item):
            grouped_files[path_item] = [path_item]
        else:
            print(f"[!] Warning: {path_item} is not a valid file or directory. Skipping.")
    
    return grouped_files


def main():
    parser = argparse.ArgumentParser(
        description="Docking Simplifier Tool - A comprehensive orchestrator for Torsion Strain analysis and Conformer generation.",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="Example of use:\n"
               "  python dock_tool.py library.mol2 -a -g -n 50\n"
               "  python dock_tool.py input.db2 --no-matrix\n"
               "  python dock_tool.py ligand.sdf -a -g\n\n"
    )
    
    parser.add_argument("input", nargs='+', help="One or more paths (files or directories)")
    
    parser.add_argument("-a", "--analyze", action="store_true", help="Execute the analysis and create the CSVs")
    parser.add_argument("-g", "--generate", action="store_true", help="Execute the generation of coordinate-updated structures")
    parser.add_argument("-n", "--num", type=int, default=-1, help="Number of Generation")
    parser.add_argument("--no-matrix", action="store_false", dest="matrix", help="Disable the creation of matrix.csv")
    parser.add_argument("-ga", "--gacha", action="store_true", help="Enable Gacha mode for random conformer generation")
    parser.add_argument("--strain", action="store_true", help="Rank generated conformers by histogram strain energy instead of MMFF94 global energy")
    parser.add_argument("--force", action="store_true", help="Ignore the internal 10k conformer generation cap")
    parser.set_defaults(matrix=True)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    
    grouped_input_files = collect_input_files(args.input)

    total_files = sum(len(files) for files in grouped_input_files.values())
    if total_files == 0:
        print("[!] ERROR: No valid supported files found.")
        sys.exit(1)
    
    # Create shared generated_conformers directory if generation is enabled
    if args.generate:
        os.makedirs("generated_conformers", exist_ok=True)
    
    with tqdm(total=total_files, desc="Processing files", position=0, leave=False, dynamic_ncols=True) as pbar:
        for input_path, file_paths in grouped_input_files.items():
            if not file_paths:
                continue

            is_directory_input = os.path.isdir(input_path)
            combined_analysis_results = [] if is_directory_input and args.analyze else None

            for file_path in file_paths:
                pbar.set_postfix_str(f"File: {os.path.basename(file_path)}")
                normalized_source_path = file_path
                temp_sdf_path = None

                try:
                    normalized_source_path, names, ms, temp_sdf_path = load_normalized_molecules(file_path)

                    for name in names:
                        m = ms.get(name)
                        if m:
                            try:
                                m.UpdatePropertyCache(strict=False)
                                Chem.FastFindRings(m)
                            except:
                                try:
                                    m.GetRingInfo().NumRings()
                                except:
                                    pass

                    positions_dict = {n: m.GetConformer().GetPositions() for n, m in ms.items()
                                      if m and m.GetNumConformers() > 0}

                    file_analysis_results = []
                    for name in names:
                        mol = ms.get(name)
                        if not mol:
                            continue

                        analysis_result = TPA.Torsion_Strain(name, mol, positions_dict)
                        if len(analysis_result) < 2 or analysis_result[1] == "NA":
                            continue

                        if args.analyze:
                            record = {
                                "name": name,
                                "analysis": analysis_result,
                                "molecule": mol,
                                "source_file": file_path,
                            }
                            if is_directory_input:
                                combined_analysis_results.append(record)
                            else:
                                file_analysis_results.append(record)

                        if args.generate:
                            input_file_base = os.path.basename(file_path)
                            generator.mol2_generator(
                                mol,
                                analysis_result,
                                args.num,
                                name,
                                args.gacha,
                                input_file_base=input_file_base,
                                strain_mode=args.strain,
                                source_file_path=normalized_source_path,
                                force=args.force,
                            )

                    if args.analyze and not is_directory_input and file_analysis_results:
                        output_dir = os.path.join(os.path.dirname(file_path), "analysis_results")
                        os.makedirs(output_dir, exist_ok=True)
                        base_name = os.path.basename(file_path)
                        new_output_path = os.path.join(output_dir, base_name)
                        TSW.generate_csv_reports(new_output_path, file_analysis_results, create_matrix=args.matrix, disable_progress=True)
                finally:
                    if temp_sdf_path and os.path.exists(temp_sdf_path):
                        try:
                            os.remove(temp_sdf_path)
                        except OSError:
                            pass

                pbar.update(1)

            if args.analyze and is_directory_input and combined_analysis_results:
                output_dir = os.path.join(input_path, "analysis_results")
                os.makedirs(output_dir, exist_ok=True)
                base_name = os.path.basename(os.path.normpath(input_path))
                new_output_path = os.path.join(output_dir, base_name + ".mol2")
                TSW.generate_csv_reports(new_output_path, combined_analysis_results, create_matrix=args.matrix, disable_progress=True)

    print("\n[*] All tasks completed. Done.")

if __name__ == "__main__":
    main()