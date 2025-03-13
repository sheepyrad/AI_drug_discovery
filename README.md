# De Novo Drug Discovery Pipeline

A comprehensive modular system for computational drug discovery using generative models, retrosynthesis, and molecular docking.

## Overview

This project integrates cutting-edge machine learning and computational chemistry to enable accelerated structure-based drug design (SBDD). The pipeline combines generative models for de novo molecule creation, retrosynthetic analysis, medicinal chemistry filtering, and automated molecular docking to identify promising drug candidates.

## Key Features

- **AI-Powered Ligand Generation**: Uses DiffSBDD, a diffusion-based generative model to create novel molecules tailored to specific protein binding sites
- **Retrosynthetic Analysis**: Implements Synformer to evaluate synthetic accessibility and generate plausible synthesis routes
- **Medicinal Chemistry Filters**: Applies drug-likeness and toxicity filters to prioritize promising compounds
- **Molecular Docking**: Supports multiple docking engines with automated pose generation and scoring
- **Multi-round Capability**: Supports iterative workflows for compound refinement across multiple rounds
- **Comprehensive Tracking**: Tracks compounds throughout the pipeline for full provenance and analysis
- **Modular Architecture**: Easily extendable with new components or custom modifications

## Pipeline Workflows

The project offers multiple workflows to accommodate different research needs:

### Standard Pipeline (`pipeline.py`)

Complete workflow with all analysis steps:
1. Ligand Generation
2. Energy Minimization (for multiple-ligand mode)
3. Pose Evaluation (concatenation + PoseBuster filtering)
4. Optional MedChem Filtering
5. Redocking

### Quick Pipeline (`pipeline_quick.py`)

Streamlined workflow for faster results:
1. Ligand Generation
2. Conversion to SMILES
3. Retrosynthesis Analysis
4. Variant Extraction
5. MedChem Filtering
6. Redocking

### Multi-round Pipeline (`pipeline_quick_multiround.py`)

Iterative exploration with multiple rounds of compound generation and analysis:
- Supports multiple sequential rounds
- Maintains compounds across rounds with unique barcoding
- Generates master reports aggregating data from all rounds

## Installation

### Prerequisites

- Linux-based operating system (Ubuntu 20.04+ recommended)
- CUDA-capable GPU (for accelerated model inference)
- Anaconda or Miniconda

### Environment Setup

1. Clone the repository:
```bash
git clone XXX
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

# Run quick pipeline (1 round)
python pipeline_quick.py --out_dir output

# Run multi-round pipeline
python pipeline_quick_multiround.py --out_dir output --num_rounds 3 --n_samples 100
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
- `--top_n`: Number of top compounds to select for analysis
- `--num_rounds`: Number of rounds for multi-round pipeline

### Example Workflow

```bash
# Generate compounds and perform full analysis
python pipeline.py --out_dir results --checkpoint models/diffsbdd_model.pt \
  --pdbfile input/target_protein.pdb --resi_list "A:28 A:29 A:30" \
  --n_samples 100 --sanitize --program_choice "qvina" \
  --center 114.817 75.602 82.416 --box_size 38 70 58
```

## Module Descriptions

- **DiffSBDD/**: Diffusion-based generative model for 3D molecule generation
- **synformer/**: Transformer-based model for retrosynthetic analysis
- **VFU/**: Virtual Flow Unity docking toolkit
- **ligand_generation.py**: Interface to DiffSBDD for molecule generation
- **retrosynformer.py**: Interface to Synformer for retrosynthetic analysis
- **medchem_filter.py**: Implementation of medicinal chemistry filters
- **energy_minimization.py**: Energy minimization and conformer generation
- **pose_evaluation.py**: Evaluation of docking poses
- **redocking.py**: Interface to docking engines
- **dock_synformer_compounds.py**: Specialized docking for Synformer compounds

## Directory Structure

The pipeline creates the following directory structure for output:

```
output_dir/
├── ligands/           # Generated ligands (SDF format)
├── minimized/         # Energy-minimized structures
├── eval/              # Pose evaluation results
├── retro/             # Retrosynthesis results
├── variants/          # Extracted retrosynthesis variants
├── filtered/          # Compounds after MedChem filtering
├── redock/            # Redocking results and poses
└── tracking/          # Tracking reports and analysis
```

## Code Style Guidelines

- **Imports**: Standard library imports first, then third-party packages, finally local modules
- **Formatting**: 4-space indentation, 120 character line limit
- **Types**: Use type hints in function signatures (Python 3.10 compatible)
- **Naming**: 
  - snake_case for functions, variables, and modules
  - CamelCase for classes
  - UPPER_CASE for constants
- **Error handling**: Use try/except blocks with specific exception types
- **Logging**: Use the Python logging module with appropriate log levels
- **Documentation**: Docstrings in triple quotes, with function parameters documented

## External Packages

**IMPORTANT**: Do not modify files in the external packages (VFU, DiffSBDD, and synformer) as they are essential external components of the pipeline, and modifications risk reproducibility.

## Testing


## Citation

If you use this pipeline in your research, please cite the following works:

[List of relevant papers and their citations]

## License

[Appropriate license information]

## Contributing

Guidelines for contributing to the project, including code style, testing requirements, and pull request process.

## Acknowledgments

This project builds upon several open-source tools and research. We gratefully acknowledge their contributions. 