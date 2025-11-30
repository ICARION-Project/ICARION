# ICARION

**Ion Collision And Reaction Instrument mOdeler for mass spectrometry and ioN mobility**

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)]()
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.0-blue)]()
[![CUDA](https://img.shields.io/badge/CUDA-12.0%2B-76B900?logo=nvidia)]()

A high-performance, GPU-accelerated framework for simulating ion trajectories in mass spectrometry and ion mobility instrumentation, with comprehensive collision physics and space charge modeling.

---

## 📋 Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Scientific Background](#scientific-background)
- [Quick Start](#quick-start)
- [Documentation](#documentation)
- [Performance](#performance)
- [Citation](#citation)
- [Contributing](#contributing)
- [License](#license)

---

## 🔬 Overview

ICARION is a modular simulation framework designed for quantitative modeling of ion transport in analytical chemistry instrumentation. It combines accurate collision physics (EHSS, Hard-Sphere, Langevin models), field-driven dynamics (electric, magnetic, RF), and computational efficiency through GPU acceleration.

**Target Applications:**
- **Ion Mobility Spectrometry (IMS)**: Drift tube, TWIMS, FAIMS
- **Mass Spectrometry**: Ion traps (Linear, 3D, Orbitrap), TOF, FT-ICR, Quadrupole filters
- **Hybrid Instruments**: IMS-MS coupling, IM-IM tandem separations
- **Method Development**: Collision cross-section (CCS) prediction, space charge optimization

**Design Philosophy:**
- ✅ **Physics-first**: Validated collision models, accurate field solvers
- ✅ **Performance**: GPU acceleration (10-100× speedup), OpenMP parallelization
- ✅ **Modularity**: Plugin architecture for forces, integrators, collision models
- ✅ **Reproducibility**: JSON configuration, HDF5 output, deterministic RNG

---

## ⚡ Key Features

### Physics Models
- **Collision Models**: 
  - EHSS (Elastic Hard-Sphere Scattering) with velocity-dependent cross-sections
  - Hard-Sphere with thermal averaging
  - Langevin polarization model
  - NoCollisions (field-only dynamics)
- **Forces**:
  - Electric fields (DC, RF, waveforms, field arrays)
  - Magnetic fields (uniform, gradient, quadrupole)
  - Space charge (O(N²) direct, O(N log N) grid, GPU P³M)
  - Hard-sphere repulsion
- **Chemistry**:
  - Ion-molecule reactions (charge transfer, proton transfer)
  - Temperature-dependent rate constants
  - Species databases (>100 ions, neutrals)

### Numerical Methods
- **Integrators**: RK4, adaptive RK45, Boris (symplectic)
- **Field Solvers**: Grid-based Poisson, field array interpolation
- **Space Charge**: FFT-based P³M algorithm (GPU), multigrid Poisson (CPU)
- **Collision Sampling**: Event-driven, null-collision methods

### GPU Acceleration (Phase 12/13) 🚀
- **Integrators**: RK4, RK45, Boris (5-50× speedup)
- **Space Charge**: P³M solver with cuFFT (10-40× speedup, 100k+ ions)
- **Thermalization**: EHSS GPU kernels (4-6× speedup)
- **Automatic Dispatch**: Threshold-based GPU/CPU selection

### Instruments (Pre-configured)
- Ion Mobility: Drift tube, traveling wave, FAIMS
- Ion Traps: Linear quadrupole (LQIT), 3D Paul, Orbitrap
- TOF: Linear, reflectron configurations
- Quadrupole: Mass filter, collision cell
- FT-ICR: Penning trap with magnetic field

---

## 📚 Scientific Background

### Collision Cross-Sections (CCS)
ICARION uses the **Trajectory Method** to compute CCS from first-principles MD simulations:

```
σ = lim(b_max→∞) πb_max² · P_scatter(b_max)
```

where `P_scatter` is the scattering probability for impact parameter `b`. Supports:
- **Velocity-dependent CCS**: EHSS model captures energy transfer
- **Rotationally-averaged CCS**: Isotropic sampling of ion orientation
- **Temperature effects**: Thermal velocity distributions (300-1000 K)

**Validation**: CCS values for calibrants (e.g., polyalanine, tetraalkylammonium) agree with literature within 2-5%.

### Space Charge Modeling
Three-level approach for computational efficiency:

| Method | Complexity | Use Case | Implementation |
|--------|-----------|----------|----------------|
| **Direct N-body** | O(N²) | N < 1k, exact | CPU, OpenMP |
| **Grid Poisson** | O(N log N) | 1k-10k ions | CPU, FFT |
| **GPU P³M** | O(N log N) | 10k-1M ions | cuFFT, CIC |

P³M achieves **333× speedup** (N=100k) vs. direct summation while maintaining ~20% accuracy for near-field interactions.

### RF Field Dynamics
Time-dependent fields (Paul traps, TWIMS) handled via:
- **Inline waveform evaluation**: `V(t) = V0 + V_RF·sin(ωt + φ)`
- **Pseudopotential approximation**: For high ω/Ω ratios
- **Adaptive timestep**: RK45 ensures numerical stability in stiff RF potentials

---

## 📦 Requirements

### Minimum Requirements
- **OS**: Linux (Ubuntu 20.04+, Debian 11+), macOS 11+, Windows 10/11 (WSL2)
- **Compiler**: GCC 9+ or Clang 10+ (C++17 support)
- **CMake**: 3.18 or higher
- **RAM**: 4 GB (CPU-only), 8 GB recommended
- **Disk**: 500 MB (source + build)

### Required Dependencies
| Library | Version | Purpose |
|---------|---------|--------|
| **HDF5** | 1.10+ | Trajectory output format |
| **JsonCpp** | 1.9+ | Configuration parsing |
| **OpenMP** | 4.5+ | CPU parallelization |
| **Catch2** | 2.13+ | Unit testing (optional) |

### GPU Acceleration (Optional)
- **CUDA Toolkit**: 12.0+ (11.8 minimum)
- **NVIDIA GPU**: Compute Capability 7.5+ (Turing, Ampere, Ada, Hopper)
  - Tested: RTX 3060, RTX 4090, RTX 5070 Ti, A100, H100
  - Minimum VRAM: 4 GB (recommended: 8+ GB for large simulations)
- **cuFFT**: Included with CUDA Toolkit

### Recommended System
```
CPU:  Intel i7/i9 or AMD Ryzen 7/9 (8+ cores)
RAM:  16-32 GB
GPU:  NVIDIA RTX 4070+ or A100 (12+ GB VRAM)
OS:   Ubuntu 22.04 LTS
```

---

## 🚀 Installation

### Ubuntu/Debian

#### 1. Install System Dependencies
```bash
# Update package list
sudo apt update

# Install build tools and libraries
sudo apt install -y \
  cmake \
  g++ \
  make \
  git \
  libhdf5-dev \
  libjsoncpp-dev \
  libopenmp-dev \
  libssl-dev

# Optional: Install CUDA Toolkit for GPU acceleration
# Method 1: From Ubuntu repositories (may be older version)
sudo apt install -y nvidia-cuda-toolkit

# Method 2: From NVIDIA (recommended, latest version)
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt update
sudo apt install -y cuda-toolkit-12-4

# Verify CUDA installation
nvcc --version
nvidia-smi
```

#### 2. Clone Repository
```bash
git clone https://github.com/chsch95/ICARION.git
cd ICARION
```

#### 3. Build
```bash
# Create build directory
mkdir build && cd build

# Configure (CPU-only)
cmake .. -DCMAKE_BUILD_TYPE=Release

# OR: Configure with GPU acceleration
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_GPU_ACCEL=ON \
  -DCMAKE_CUDA_ARCHITECTURES="75;80;86;89;90"  # Adjust for your GPU

# Build (parallel, ~5 minutes on 8 cores)
make -j$(nproc)

# Run tests (optional but recommended)
ctest --output-on-failure
```

#### 4. Verify Installation
```bash
# Check binary
./src/icarion_main --version

# Run quick test (should complete in ~10 seconds)
./src/icarion_main ../examples/ims_basic.json

# Check output
ls -lh ../results/ims_basic/trajectories.h5
```

### macOS

```bash
# Install Homebrew (if not already installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install cmake gcc hdf5 jsoncpp libomp openssl

# Clone and build
git clone https://github.com/chsch95/ICARION.git
cd ICARION
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(sysctl -n hw.ncpu)

# Note: GPU acceleration not available on macOS (Metal support planned for v1.1)
```

### Windows (WSL2)

```powershell
# 1. Install WSL2 and Ubuntu 22.04
wsl --install -d Ubuntu-22.04

# 2. Inside WSL2, follow Ubuntu instructions above

# 3. For GPU support, install NVIDIA CUDA on WSL
# See: https://docs.nvidia.com/cuda/wsl-user-guide/index.html
```

### Docker (Cross-platform)

```bash
# CPU-only image
docker pull ghcr.io/chsch95/icarion:1.0.0
docker run -v $(pwd)/examples:/data ghcr.io/chsch95/icarion:1.0.0 \
  icarion_main /data/ims_basic.json

# GPU-enabled image (requires nvidia-docker)
docker pull ghcr.io/chsch95/icarion:1.0.0-cuda
docker run --gpus all -v $(pwd)/examples:/data \
  ghcr.io/chsch95/icarion:1.0.0-cuda \
  icarion_main /data/ims_basic.json
```

---

## 🚀 Quick Start

### Basic Usage

### Build
```bash
git clone https://github.com/chsch95/ICARION.git
cd ICARION
mkdir build && cd build

# CPU-only build
cmake .. -DCMAKE_BUILD_TYPE=Release

# GPU-accelerated build (requires CUDA 12.0+)
cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_GPU_ACCEL=ON

make -j$(nproc)
```

### Run Example
```bash
# IMS simulation (1000 ions, EHSS collisions)
./src/icarion_main ../examples/ims_basic.json

# Output: results/ims_basic/trajectories.h5
# Analyze with Python/HDFView
```

### Configuration
Minimal JSON example:
```json
{
  "simulation": {
    "total_time_s": 1e-3,
    "dt_s": 1e-9,
    "integrator": "RK45"
  },
  "physics": {
    "collision_model": "EHSS",
    "enable_space_charge": true
  },
  "ions": {
    "species": [
      {"id": "H3O+", "count": 1000, "position": {"type": "gaussian"}}
    ]
  },
  "domains": [
    {
      "name": "drift_tube",
      "instrument": "IMS",
      "geometry": {"length_m": 0.1, "radius_m": 0.01},
      "env": {"pressure_Pa": 300, "gas_species": "N2"},
      "fields": {"DC": {"EN_Td": 10.0}}
    }
  ]
}
```

See [`examples/`](examples/) for 15+ instrument configurations.

---

## 📖 Documentation

| Document | Description |
|----------|-------------|
| [**Architecture Guide**](docs/ARCHITECTURE.md) | System design, module structure, GPU pipeline |
| [**Configuration Guide**](docs/CONFIG_GUIDE.md) | JSON schema, parameter reference |
| [**CLI Usage**](docs/CLI_USAGE.md) | Command-line options, batch processing |
| [**Collision Models**](docs/COLLISION_MODELS.md) | EHSS, Hard-Sphere, Langevin physics |
| [**HDF5 Output**](docs/HDF5_OUTPUT_STRUCTURE.md) | Data format, trajectory analysis |
| [**Developer's Guide**](docs/DEVELOPERS_GUIDE.md) | Code structure, contribution workflow |

**Key Sections:**
- [Module Structure](docs/ARCHITECTURE.md#module-structure) - Code organization
- [Force System](docs/ARCHITECTURE.md#force-system-architecture) - Plugin interface
- [GPU Acceleration](docs/ARCHITECTURE.md#gpu-acceleration-architecture) - Phase 12/13 details
- [Validation](validation/README.md) - CCS benchmarks, thermalization tests

---

## ⚡ Performance

### Benchmarks (NVIDIA RTX 5070 Ti, 16GB VRAM)

**GPU vs CPU Speedup (100k ions, 10ns simulation):**
```
Integrator (RK4):        12× faster    (60 vs 720 ms/timestep)
Space Charge (P³M):      40× faster    (50 vs 2000 ms/timestep)
EHSS Thermalization:     5× faster     (200 vs 1000 ms/10k collisions)
Overall (field-only):    8-15× faster
Overall (space charge):  20-50× faster
```

**Scaling (Space Charge P³M on GPU):**
| Ions | Grid | CPU Time | GPU Time | Speedup |
|------|------|----------|----------|---------|
| 1k   | 64³  | 20 ms    | 2 ms     | 10×     |
| 10k  | 128³ | 200 ms   | 15 ms    | 13×     |
| 100k | 256³ | 2 sec    | 50 ms    | 40×     |
| 1M   | 256³ | 3 min    | 200 ms   | 900×    |

**Memory Requirements:**
- CPU: ~200 MB (10k ions) + field grids
- GPU: ~500 MB (10k ions) + 2GB (P³M solver, 256³ grid)
- Peak: ~8 GB VRAM for 1M ions with space charge

### Parallelization
- **OpenMP**: Force computation, collision handling (linear scaling to 32 cores)
- **CUDA**: Integrators, space charge, thermalization (tested up to RTX 5090)
- **Hybrid**: CPU collision detection + GPU force evaluation

---

## 📄 Citation

If you use ICARION in your research, please cite:

```bibtex
@software{icarion2025,
  author       = {Schmidt, Christian},
  title        = {{ICARION: Ion Collision And Reaction Instrument mOdeler}},
  year         = {2025},
  publisher    = {GitHub},
  version      = {1.0.0},
  url          = {https://github.com/chsch95/ICARION}
}
```

**Related Publications:**
- *Coming soon* - CCS validation paper (J. Am. Soc. Mass Spectrom.)
- *Coming soon* - GPU acceleration methods (Int. J. Mass Spectrom.)

---

## 🤝 Contributing

Contributions welcome! See [DEVELOPERS_GUIDE.md](docs/DEVELOPERS_GUIDE.md) for:
- Code style guidelines (clang-format, naming conventions)
- Testing requirements (Catch2, >80% coverage)
- Pull request workflow
- Feature roadmap (v1.1: BEM field solver, multi-GPU)

**Areas for Contribution:**
- 🔬 **Physics**: Reaction mechanisms, CID models, ion heating
- 🚀 **Performance**: Multi-GPU, CPU vectorization, memory optimization
- 🔧 **Instruments**: SLIM, SOFI-IM, cyclic IMS configurations
- 📊 **Analysis**: Python tools, visualization, data pipelines

---

## 📜 License

ICARION is licensed under the **MIT License**.

```
MIT License

Copyright (c) 2025 ICARION Project Contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

See [LICENSE](LICENSE) for full text.

---

## 🔗 Resources

- **Homepage**: [https://github.com/chsch95/ICARION](https://github.com/chsch95/ICARION)
- **Issue Tracker**: [GitHub Issues](https://github.com/chsch95/ICARION/issues)
- **Discussions**: [GitHub Discussions](https://github.com/chsch95/ICARION/discussions)
- **Examples**: [`examples/`](examples/) - 15+ instrument configs
- **Validation**: [`validation/`](validation/) - CCS benchmarks, convergence tests

**External Links:**
- [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit) - GPU acceleration
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) - Output format
- [Catch2](https://github.com/catchorg/Catch2) - Testing framework

---

<div align="center">

**Built with ❤️ for the analytical chemistry community**

[⬆ Back to Top](#icarion)

</div>