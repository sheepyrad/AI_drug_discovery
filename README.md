# De Novo Drug Discovery Pipeline

A comprehensive modular system for computational drug discovery using generative models, retrosynthesis, and molecular docking.

## Overview

This project integrates cutting-edge machine learning and computational chemistry to enable accelerated structure-based drug design (SBDD). The pipeline combines generative models for de novo molecule creation, retrosynthetic analysis, medicinal chemistry filtering, and automated molecular docking to identify promising drug candidates.

## Project Structure

```
.
├── src/                    # External dependencies
│   ├── DiffSBDD/          # Diffusion-based generative model
│   ├── Synformer/         # Retrosynthesis analysis
│   └── VFU/               # Virtual Flow Unity for docking
├── utils/                  # Utility functions
│   ├── __init__.py        # Package initialization
│   ├── ligand_generation.py
│   ├── retrosynformer.py
│   ├── medchem_filter.py
│   ├── energy_minimization.py
│   ├── dock_synformer_compounds.py
│   ├── pose_evaluation.py
│   └── redocking.py
├── pipeline.py            # Standard pipeline script
└── pipeline_quick_multiround.py  # Multi-round quick pipeline
```

## Key Features

- **AI-Powered Ligand Generation**: Uses DiffSBDD, a diffusion-based generative model to create novel molecules tailored to specific protein binding sites
- **Retrosynthetic Analysis**: Implements Synformer to evaluate synthetic accessibility and generate plausible synthesis routes
- **Medicinal Chemistry Filters**: Applies drug-likeness and toxicity filters to prioritize promising compounds
- **Molecular Docking**: Supports multiple docking engines with automated pose generation and scoring
- **Multi-round Capability**: Supports iterative workflows for compound refinement across multiple rounds
- **Real-time Tracking**: Monitors compound progression through the pipeline with detailed status updates
- **Comprehensive Tracking**: Tracks compounds throughout the pipeline with timestamps and status changes:
  - Generation
  - Retrosynthesis
  - MedChem Filtering
  - Docking
- **Modular Architecture**: Easily extendable with new components or custom modifications

## Pipeline Workflows

### Standard Pipeline (`pipeline.py`)
Complete workflow with all analysis steps:
1. Ligand Generation
2. Energy Minimization (for multiple-ligand mode)
3. Pose Evaluation (concatenation + PoseBuster filtering)
4. Optional MedChem Filtering
5. Redocking

### Quick Multi-round Pipeline (`pipeline_quick_multiround.py`)
Streamlined workflow for iterative exploration:
1. Ligand Generation
2. Convert to SMILES
3. Retrosynthesis Analysis
4. Variant Extraction
5. MedChem Filtering
6. Redocking
7. Real-time Progress Tracking

## Installation

### Prerequisites
- Linux-based operating system (Ubuntu 20.04+ recommended)
- CUDA-capable GPU (for accelerated model inference)
- Anaconda or Miniconda

### Environment Setup
1. Clone the repository:
```bash
git clone <repository_url>
cd denovo-drug-discovery
```

2. Create and activate the conda environment:
```bash
conda env create -f environment.yml
conda activate sbdd-env-exp
```

## Usage

### Basic Usage
```bash
# Run standard pipeline
python pipeline.py --out_dir output

# Run multi-round pipeline
python pipeline_quick_multiround.py --out_dir output --num_rounds 3
```

### Key Parameters
- `--out_dir`: Directory for all output files (required)
- `--checkpoint`: Path to DiffSBDD checkpoint
- `--pdbfile`: Input protein structure file
- `--resi_list`: List of residues defining the binding site
- `--n_samples`: Number of ligands to generate
- `--sanitize`: Apply RDKit sanitization to generated molecules
- `--protein_file`: Prepared protein file for docking
- `--receptor`: Path to receptor file
- `--program_choice`: Docking program to use (default: "qvina")
- `--scoring_function`: Scoring function (default: "nnscore2")
- `--center`: Center coordinates for docking box
- `--box_size`: Size of docking box
- `--exhaustiveness`: Docking exhaustiveness
- `--num_rounds`: Number of rounds for multi-round pipeline
- `--max_variants`: Maximum variants per compound for retrosynthesis

### Output Structure
```
output/
├── master_tracking/              # Master tracking across all rounds
│   └── master_compound_tracking_report.csv
├── round_1/                      # Round-specific directories
│   ├── tracking_report.csv       # Round-specific tracking
│   ├── ligand_generation/       
│   ├── retrosyn_results/
│   ├── filter_results/
│   └── docking_results/
└── round_N/                      # Subsequent rounds...
```

## Tracking Reports
The pipeline generates detailed tracking reports that include:
- Compound IDs and barcodes
- Generation timestamps
- Processing status
- Source information
- Docking scores and poses
- Parent-child relationships for variants

## Contributing
Contributions are welcome! Please read our contributing guidelines and code of conduct.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments
- DiffSBDD team for the generative model
- Synformer team for retrosynthesis capabilities
- Virtual Flow Unity (VFU) team for docking infrastructure 