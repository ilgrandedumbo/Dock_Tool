# Docking Simplifier Tool (`dock_tool`)

A specialized, comprehensive Python orchestrator designed for molecular torsion strain analysis, automated conformer structure generation, and advanced SMARTS pattern mapping. This tool integrates Torsion Library rules with molecular topologies using RDKit to automate the evaluation of torsion strain energy parameters and geometric variations.

## Key Capabilities & Module Overview

- **Main Orchestrator (`dock_tool.py`):** The primary command-line interface that automates end-to-end processing. It seamlessly handles single molecular files or entire directories, manages file extraction (including transparent support for `.gz` compressed formats), and coordinates analysis and conformation generation pipelines.
- **Core Torsion Evaluation (`TL_Functions_modify_2.py`):** Parses specific and general hierarchy classes directly from the Torsion Library XML database. It maps exact and approximate rules onto RDKit molecular structures and implements custom central-bond resolution logic to eliminate duplicate or redundant structural matching.
- **Parallelized Batch Processing (`torsion_pattern_analyzer.py`):** Implements a high-performance multiprocessing framework (`ProcessPoolExecutor`) to calculate dihedral angles and interpolate continuous data points across thousands of chemical compounds efficiently.
- **Conformer Generation Engine (`generator_tool.py`):** Evaluates allowed combinations of optimal dihedrals and synthesizes new valid 3D conformers based on topological rules.
- **Statistical Report Generator (`Torsion_Strain_Writer.py`):** Aggregates cross-run data and generates comprehensive outputs, exporting clean summaries, raw statistics, and dihedral correlation matrices directly into tabular CSV format.

## Repository Structure

```text
dock_tool_repo/
├── LICENSE                    # Creative Commons BY-NC-SA 4.0 license file
├── README.md                  # Project documentation
├── requirements.txt           # Python software dependencies
├── dock_tool.py               # Main pipeline orchestrator script
├── TL_Functions_modify_2.py   # Torsion Library parsing and supplier utilities
├── torsion_pattern_analyzer.py# Multiprocessing strain analysis core
├── Torsion_Strain_Writer.py   # CSV report writer and matrix builder
├── generator_tool.py          # Geometry and conformer generator engine
└── TL_2.1_VERSION_6.xml       # Reference Torsion Library XML database
```
## Installation

### Prerequisites

- **Python 3.9+** is required to support the underlying typing and core library implementations.
- Using a virtual environment tool (such as conda or venv) is highly recommended to isolate environment variables.

### Setup

1. Clone this repository to your local system: ```bash git clone [https://github.com/](https://github.com/)[YourGitHubUsername]/dock_tool.git
cd dock_tool```

2. Install all required packages using the provided **requirements.txt**: ```bash pip install -r requirements.txt```

Note: The environment requires rdkit for handling chemoinformatics operations, numpy for matrix calculations, and tqdm for terminal progress visualization.

## Usage

Run the main pipeline by feeding the orchestrator a molecular structure file (.smi, .sdf, .mol2, or .db2). The interface automatically detects extension types and controls execution flow:  ```bash python dock_tool.py /path/to/input_file.mol2 [options] ```

### Command-Line Arguments

For all optional argumets types ```bash python dock_tool.py -h```


## License & Terms of Use

This project is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA 4.0).

- Academic & Research Use: Completely free for students, academic researchers, universities, and non-profit institutions.

- Commercial & SaaS Restriction: Commercial exploitation, integration into proprietary closed-source applications, or deployment as a paid web service or Cloud SaaS platform is strictly prohibited without explicit, prior written authorization from the authors.

### Third-Party Components

This repository includes molecule parsing utilities and DB2 chemical format suppliers adapted from DOCK3.7 (Copyright 2025 The Regents of the University of California), originally released under the permissive BSD-3-Clause license. The original copyright notices and full liability disclaimers are preserved intact inside the source header of TL_Functions_modify_2.py.

## Authors & Citation

Developed by Niccolò Papini and Sardelli.

If you use this orchestrator or its core analytical modules in a scientific publication, please cite this repository:

```text 
Papini, N.; Sardelli. Docking Simplifier Tool (dock_tool), 2026. GitHub Repository.```
