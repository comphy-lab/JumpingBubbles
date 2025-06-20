# Jumping Bubbles

[![DOI](https://zenodo.org/badge/744202007.svg)](https://doi.org/10.5281/zenodo.14602622)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![OpenMP](https://img.shields.io/badge/OpenMP-enabled-brightgreen.svg)](https://www.openmp.org/)
[![MPI](https://img.shields.io/badge/MPI-enabled-brightgreen.svg)](https://www.open-mpi.org/)
[![GitHub issues](https://img.shields.io/github/issues/VatsalSy/jumping-bubbles)](https://github.com/VatsalSy/jumping-bubbles/issues)
[![GitHub pull requests](https://img.shields.io/github/issues-pr/VatsalSy/jumping-bubbles)](https://github.com/VatsalSy/jumping-bubbles/pulls)
[![GitHub release](https://img.shields.io/github/v/release/VatsalSy/jumping-bubbles)](https://github.com/VatsalSy/jumping-bubbles/releases)
[![Documentation Status](https://img.shields.io/badge/docs-latest-brightgreen.svg)](http://basilisk.fr/sandbox/README)

A high-performance computational framework for studying bubble coalescence and jumping phenomena on substrates using [Basilisk C](http://basilisk.fr/). This repository provides a complete suite of simulation tools, post-processing utilities, and validated test cases for investigating two-phase free-surface flows with adaptive mesh refinement.

## Key Features

### **Advanced Two-Phase Flow Modeling**
- Volume of Fluid (VOF) method for sharp interface tracking
- Conservative momentum advection coupled with VOF field
- Accurate surface tension implementation using the Brackbill method
- Reduced gravity approach for specific physical scenarios. This keeps the system of equations well-balanced.

### **High-Performance Computing**
- Adaptive octree mesh refinement based on multiple criteria:
  - Interface fraction gradients
  - Local curvature
  - Velocity field variations
- Parallel computation support:
  - OpenMP for shared-memory systems
  - MPI for distributed computing (HPC clusters)

### **Comprehensive Analysis Tools**
- Real-time monitoring of physical quantities
- Advanced post-processing capabilities:
  - 2D/3D visualization scripts
  - Interface reconstruction
  - Energy and momentum analysis
  - Trajectory tracking

## Project Structure

```
├── basilisk/src/               # Core Basilisk C framework
│   ├── navier-stokes/         # Flow solvers
│   │   ├── centered.h         # Main centered NS solver
│   │   └── conserving.h       # Conservative form solver
│   ├── two-phase*.h           # Two-phase flow implementations
│   ├── vof.h                  # Volume of Fluid method
│   ├── curvature.h           # Interface property calculations
│   ├── tension.h             # Surface tension (Brackbill)
│   ├── integral.h            # Surface tension integral form
│   ├── viscosity.h           # Implicit viscous stress solver
│   └── axi.h                 # Axisymmetric computations
├── src-local/                 # Project-specific modifications
├── postProcess/               # Analysis and visualization
│   ├── Video2DSlice.py       # 2D slice animations
│   ├── Video3D.py            # 3D visualization
│   ├── Visualization3D.ipynb  # Interactive 3D notebook
│   └── getFacets3D.c         # Interface extraction
└── simulationCases/                # Validation cases
    ├── JumpingBubbles.c      # Standard test case
    ├── JumpingBubbles-hydrophilic.c  # Hydrophilic substrate
    └── runCases*.sh          # Execution scripts
```

## Installation

### Prerequisites
- GCC compiler (version 7.0+)
- MPI implementation (OpenMPI or MPICH)
- Python 3.7+ (for post-processing)
- Required Python packages: numpy, matplotlib, scipy
- Xcode Command Line Tools (for MacOS users)

### Setting up Basilisk
1. Clone this repository:
   ```bash
   git clone https://github.com/VatsalSy/jumping-bubbles.git
   cd jumping-bubbles
   ```

2. Run the installation script:
   ```bash
   ./reset_install_requirements.sh
   ```
   This script will:
   - Install Basilisk if not present (use `--hard` flag for fresh installation)
   - Configure the environment automatically
   - Create `.project_config` with required paths
   - Verify the installation

3. Verify installation:
   ```bash
   source .project_config
   qcc --version
   ```

## Usage

### Running Simulations

#### Local Development (OpenMP)
```bash
# Compile
qcc -O2 -Wall -disable-dimensions -fopenmp -I$(PWD)/src-local testCases/JumpingBubbles.c -o JumpingBubbles -lm

# Run with 4 threads
export OMP_NUM_THREADS=4
./JumpingBubbles
```

#### HPC Deployment (MPI)
1. Generate initial condition locally first
2. Compile with MPI support:
   ```bash
   CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions \
   -I$(PWD)/src-local testCases/JumpingBubbles.c -o JumpingBubbles -lm
   ```
3. Use provided Slurm script: `testCases/runSnellius.sbatch`

### Post-Processing

#### Visualization
1. 2D Slice Videos:
   ```bash
   python postProcess/Video2DSlice.py <simulation_directory>
   ```

2. 3D Visualization:
   ```bash
   python postProcess/Video3D.py <simulation_directory>
   ```

3. Interactive Analysis:
   - Open `postProcess/Visualization3D.ipynb` in Jupyter

#### Data Analysis
- Use `getFacets3D.c` to extract interface geometry
- Various slice extraction tools available in `postProcess/`

## Contributing

We welcome contributions to improve Jumping Bubbles! Here's how you can help:

### Reporting Issues
Click on the "Issues" tab above or use these quick links:
- [Report a Bug](https://github.com/VatsalSy/jumping-bubbles/issues/new?template=bug_report.md&labels=bug)
- [Request a Feature](https://github.com/VatsalSy/jumping-bubbles/issues/new?template=feature_request.md&labels=enhancement)
- [Ask a Question](https://github.com/VatsalSy/jumping-bubbles/issues/new?template=question.md&labels=question)
- [Open a Blank Issue](https://github.com/VatsalSy/jumping-bubbles/issues/new)

### Making Changes
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Make your changes
4. Run tests and ensure they pass
5. Update documentation if needed
6. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
7. Push to the branch (`git push origin feature/AmazingFeature`)
8. Open a Pull Request

Please ensure your PR:
- Clearly describes the changes
- Includes any relevant tests
- Updates documentation as needed
- References any related issues

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Citation
If you use this code in your research, please cite:
```bibtex
@software{jumping_bubbles_2024,
  author       = {Sanjay, V. and Yang, R.},
  title        = {Jumping Bubbles: A Computational Framework for Studying Bubble Coalescence},
  year         = {2025},
  publisher    = {Zenodo},
  version      = {v2.0},
  doi          = {10.5281/zenodo.14602622},
  url          = {https://doi.org/10.5281/zenodo.14602622}
}
```

## Acknowledgments
- Based on the [Basilisk C](http://basilisk.fr/) framework.