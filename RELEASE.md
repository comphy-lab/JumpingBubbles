# Release Notes

## Version 2.0
[![Version](https://img.shields.io/badge/version-v2.0-brightgreen.svg)](https://github.com/comphy-lab/JumpingBubbles/releases)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/744202007.svg)](https://doi.org/10.5281/zenodo.14602622)
[![OpenMP](https://img.shields.io/badge/OpenMP-enabled-brightgreen.svg)](https://www.openmp.org/)
[![MPI](https://img.shields.io/badge/MPI-enabled-brightgreen.svg)](https://www.open-mpi.org/)
[![Python](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

**Release Date:** June 20, 2025

### 🚀 Major Highlights

This release represents a significant evolution of the Jumping Bubbles framework, featuring a complete documentation overhaul, improved code organization, enhanced visualization capabilities, and numerous quality-of-life improvements.

### 📚 Documentation Overhaul

- **New Documentation Website**
  - Comprehensive GitHub Pages documentation with modern, responsive design
  - Interactive search functionality with command palette
  - Dark/light theme toggle support
  - Mobile-friendly interface with custom styling
  - Full HTML documentation for all source files

- **Enhanced Content**
  - Detailed inline documentation added to all source files
  - Comprehensive installation and usage guides
  - Added contribution guidelines and issue templates
  - Improved code examples and tutorials

### 🏗️ Project Restructuring

- **Reorganized Directory Structure**
  - Simulation files moved to `simulationCases/` directory
  - Post-processing tools consolidated in `postProcess/` folder
  - Clear separation between core simulation code and utilities
  - Better naming conventions throughout the project

### ✨ New Features

- **Enhanced Contact Angle Support**
  - Extended implementation for arbitrary contact angles
  - Improved boundary condition handling for contact lines
  - Added Bond number variable for better physical modeling

- **Hydrophilic Surface Modeling**
  - New `JumpingBubbles-hydrophilic.c` variant
  - Enables comparative studies between different wetting conditions

- **Advanced Visualization**
  - Updated `Visualization3D.ipynb` Jupyter notebook
  - Enhanced 2D/3D video generation capabilities
  - PNG image support for documentation

### 🐛 Bug Fixes

- Fixed boolean parsing in data extraction utilities (`getDataXSlice.c`, `getDataZSlice.c`)
- Resolved NameError in `Video3D.py` subprocess exception handling
- Fixed subprocess import alias in `Video2DSlice.py`
- Corrected unescaped dash in regex patterns
- Fixed YAML indentation in GitHub workflows
- Added localStorage error handling for theme toggle
- Added guard checks for window.searchHelper

### 🔧 Improvements

- **Code Quality**
  - Standardized naming conventions across postProcess folder
  - Improved error handling in Python visualization scripts
  - Enhanced .gitignore patterns for better repository management
  - Consistent code formatting and documentation style

- **Build System**
  - GitHub Actions workflows for automated documentation deployment
  - Search index generation workflow
  - Improved installation script with better Basilisk handling

- **HPC Support**
  - Updated SLURM batch scripts
  - Better handling of distributed computing scenarios
  - Maintained OpenMP/MPI compatibility

### 📊 Known Issues

- Basilisk's `distance.h` remains incompatible with MPI (unchanged from v1.0)
- Some HPC environments still require manual compiler configuration

### 🔮 Future Directions

- Automated regression testing framework implementation
- Further optimization of parallel computing performance
- Enhanced interface tracking algorithms
- Additional test cases for validation

### 📝 Citation

```bibtex
@software{jumping_bubbles_2025,
  author       = {Sanjay, V. and Yang, R.},
  title        = {Jumping Bubbles: A Computational Framework for Studying Bubble Coalescence},
  year         = {2025},
  publisher    = {Zenodo},
  version      = {v2.0},
  doi          = {10.5281/zenodo.14602622},
  url          = {https://doi.org/10.5281/zenodo.14602622}
}
```

---

## Version 1.0
[![Version](https://img.shields.io/badge/version-v1.0-brightgreen.svg)](https://github.com/comphy-lab/JumpingBubbles/releases)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/744202007.svg)](https://doi.org/10.5281/zenodo.14602622)
[![OpenMP](https://img.shields.io/badge/OpenMP-enabled-brightgreen.svg)](https://www.openmp.org/)
[![MPI](https://img.shields.io/badge/MPI-enabled-brightgreen.svg)](https://www.open-mpi.org/)
[![Python](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

**Release Date:** January 5, 2025

### Highlights

#### Expanded Two-Phase Flow Capabilities
- VOF-based interface tracking with an option for the conservative momentum advection module
- Added the reduced gravity approach for precise control over effective body forces

#### Improved HPC Support
- MPI integration for large-scale parallel runs on HPC clusters
- Refined instructions and Slurm scripts for streamlined deployment on typical HPC infrastructures

#### Enhanced Post-Processing Tools
- New Python scripts for 2D and 3D visualization (`Video2DSlice.py`, `Video3D.py`)
- Jupyter notebook (`Visualization3D.ipynb`) for interactive exploration of simulation data
- `getFacets3D.c` for interface extraction and geometric analysis

#### Robust Initialization Workflow
- Option to generate initial conditions in a single-threaded (OpenMP) run before migrating to MPI
- Automatic environment configuration via `.project_config` to minimize setup inconsistencies

### Bug Fixes
- Resolved mesh refinement edge cases when using `distance.h`
- Corrected velocity boundary conditions on domain walls for axisymmetric setups
- Addressed race conditions in OpenMP runs by clarifying memory allocation routines

### Known Issues
- Basilisk's `distance.h` remains incompatible with MPI, necessitating local initialization before parallel runs
- Some HPC environments require manual customization of compiler and library paths in `.project_config`

### Future Directions
- Addition of automated regression tests in `testCases/` for routine CI checks
- Enhancements to post-processing scripts for slice-by-slice interface tracking
- Further optimization of load balancing in MPI runs

### Citation

```bibtex
@software{jumping_bubbles_2024,
  author       = {Sanjay, V. and Yang, R.},
  title        = {Jumping Bubbles: A Computational Framework for Studying Bubble Coalescence},
  year         = {2024},
  publisher    = {Zenodo},
  version      = {v1.0},
  doi          = {10.5281/zenodo.14602622},
  url          = {https://doi.org/10.5281/zenodo.14602622}
}
```