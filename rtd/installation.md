# Installation

ICARION can be used either from prebuilt release packages or by building the source code.

## Prebuilt releases

For most users, the easiest starting point is the latest release package from the GitHub repository:

<https://github.com/ICARION-Project/ICARION/releases/latest>

Download the package for your operating system and follow the release notes and the 
`ICARION-Launcher-Guide.md`. The package layouts and launcher 
entry points are summarized in
[Release packages and launcher](release-packages.md).

## Build from source

A source build is recommended for developers, for users who want to inspect or modify the code, and for systems where no prebuilt package is available.

### Requirements

A typical CPU build requires:

- a C++17 compiler,
- CMake,
- HDF5,
- Eigen3,
- jsoncpp,
- nlohmann_json,
- OpenSSL,
- BLAS,
- spdlog,
- cxxopts,
- and optionally OpenMP.

On Ubuntu or WSL, the basic dependencies can be installed with:

```bash
sudo apt update
sudo apt install -y \
  build-essential cmake git pkg-config \
  libeigen3-dev libjsoncpp-dev nlohmann-json3-dev libhdf5-dev \
  libssl-dev libopenblas-dev libspdlog-dev libcxxopts-dev
```

Then clone and build ICARION:

```bash
git clone https://github.com/ICARION-Project/ICARION.git
cd ICARION

cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
```

Run a basic example:

```bash
./build/src/icarion_main examples/ims/ims_basic.json
```

## Optional GPU build

GPU-related code can be compiled if CUDA is available:

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DUSE_GPU_ACCEL=ON
cmake --build build -j"$(nproc)"
```

!!! warning
    GPU support should be treated as experimental unless the release notes for the version you are using explicitly state otherwise.

After installation, continue with [First simulation](first-simulation.md). For
command-line flags, see [CLI reference](cli-reference.md). To verify a build, see
[Validation](validation.md).
