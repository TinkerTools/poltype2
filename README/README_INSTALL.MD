## Installation

### Prerequisites

| Program | Version | Notes |
|---------|---------|-------|
| Python | >= 3.9 | With pip |
| Tinker | >= 8.9.4 | CPU build required |
| GDMA | >= 2.3 | Distributed multipole analysis |
| Psi4 **or** Gaussian 09/16 | — | Quantum chemistry engine |

### 1. Install Poltype2

```bash
git clone https://github.com/TinkerTools/poltype2.git
cd poltype2
pip install -e ".[dev]"
```

This installs `poltype2` as an editable package with the `poltype2` CLI
command.  Core dependencies (RDKit, NumPy) are pulled in automatically.

### 2. Install Tinker (CPU)

```bash
git clone git@github.com:TinkerTools/Tinker.git tinker
cd tinker/fftw
./configure --prefix=$(pwd) --enable-openmp --enable-threads
make -j 8 && make install
cd ..
cp ./make/Makefile ./source
cd source
```

Open `Makefile`, comment out the default OS section (Mac) and uncomment the
Linux lines.  Keep the `RENAME` line commented out.  Then:

```bash
make -j 8
```

Add the resulting executables (`analyze`, `poledit`, `potential`, `minimize`,
`dynamic`, `bar`) to your `PATH`.

#### Optional: Tinker-GPU

Tinker-GPU is **not** required for parameterisation.  See the
[Tinker9 repository](https://github.com/TinkerTools/tinker9) for GPU build
instructions.

### 3. Install a Quantum Chemistry Package

**Psi4 (recommended, free licence):**

```bash
conda install -c conda-forge psi4
```

**Gaussian 09/16:** Follow your site's installation instructions and ensure
`g09` or `g16` is on your `PATH`.

**PySCF (lightweight alternative):**

```bash
pip install pyscf
```

### 4. Install GDMA

Download GDMA >= 2.3 from the official source and add the `gdma` binary to
your `PATH`.

### 5. Optional: Install FeNNol (ML Potentials)

```bash
git clone https://github.com/thomasple/FeNNol.git
cd FeNNol && pip install .
```

### 6. Environment Variables

Before running Poltype, ensure these programs are discoverable:

```bash
# Example: source this before running poltype
export PATH=/path/to/tinker/bin:$PATH
export PATH=/path/to/gdma/bin:$PATH
export PSI_SCRATCH=/scratch/$USER/
# If using Gaussian:
export g09root=/opt/gaussian
source $g09root/g09/bsd/g09.profile
export GAUSS_SCRDIR=/scratch/$USER/
```

> **Tip:** Do not put `source` commands for multiple conflicting environments
> in the same `.bashrc` file.

### 7. Verify Installation

```bash
poltype2 --version
# or
python -m poltype --version
```
