# Changelog

## [2024-03-13] Real-time Tracking Implementation

### Added
- New `update_tracking_report` function for real-time compound tracking
- Real-time status updates for compounds through their lifecycle:
  - "GENERATED": Initial compound creation
  - "SYNTHETIZED": After retrosynthesis
  - "FILTERED_PASS": After passing MedChem filters
  - "DOCKED": After successful docking
- Timestamps for all tracking entries
- Round-specific tracking reports
- Master tracking report across all rounds

### Modified
- Refactored pipeline_quick_multiround.py for real-time tracking:
  - Removed batch processing in favor of immediate updates
  - Integrated tracking updates at each step of compound processing
  - Streamlined compound processing flow
  - Added better error handling and logging
- Enhanced logging to provide more detailed progress information

### Improvements
- Real-time visibility into pipeline progress
- Lower memory usage (no data accumulation)
- Better resilience to failures (immediate data persistence)
- Enhanced compound lifecycle tracking
- Improved monitoring capabilities for long-running pipelines

## [2024-03-13] Project Restructuring

### Added
- Created `src/` directory for external dependencies
- Created `utils/` directory for utility scripts
- Added this CHANGELOG.md file
- Created `utils/__init__.py` to make utils a proper Python package

### Moved
- External packages moved to `src/`:
  - DiffSBDD
  - Synformer
  - VFU

- Utility scripts moved to `utils/`:
  - ligand_generation.py
  - retrosynformer.py
  - medchem_filter.py
  - energy_minimization.py
  - dock_synformer_compounds.py
  - energy_minimization_module.py
  - pose_evaluation.py
  - redocking.py

### Modified
- Updated import statements in pipeline.py to use utils package
- Updated import statements in pipeline_quick_multiround.py to use utils package
- Updated paths in utility scripts to reference external packages in src/:
  - ligand_generation.py: Updated DiffSBDD path
  - retrosynformer.py: Updated Synformer paths
  - redocking.py: Updated VFU paths and fixed input directory path

### Pipeline Scripts (Unchanged Location)
The following scripts remain in the root directory:
- pipeline.py
- pipeline_quick_multiround.py

### Notes
- This restructuring improves project organization by:
  - Separating external dependencies from internal code
  - Grouping utility functions in a dedicated directory
  - Keeping main pipeline scripts easily accessible in the root directory
  - Making the utils directory a proper Python package
- All utility scripts have been updated to properly reference external packages in their new location under src/ 