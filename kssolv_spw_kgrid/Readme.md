# KSSOLV-SPW

KSSOLV-SPW is an extension of KSSOLV that implements the symmetrized
plane-wave (SPW) method for accelerating iterative diagonalization in
plane-wave density functional theory (DFT) using space-group symmetry.

The implementation retains the conventional KSSOLV self-consistent-field
(SCF) framework and block-Davidson eigensolver, while decomposing the
plane-wave eigenproblem into independent symmetry-adapted subproblems.

This repository contains the implementation used in the associated
article:

> **Symmetrized plane waves for accelerating iterative diagonalization
> in density functional theory**

## Repository structure

The example calculations are located in

```text
kssolv_spw_kgrid/
├── ...
└── spw_example/
    ├── Gamma/
    ├── Kpoint/
    ├── structure_spg/
    └── Reproducing the SPW benchmark examples.md
```

### `spw_example/Gamma`

Contains the scripts for reproducing the representative
$\Gamma$-point calculations reported in the article.

The main test script is

```text
spw_example/Gamma/test_spw.m
```

Running `test_spw` performs, in sequence:

1. a conventional SCF calculation using the target number of bands;
2. a conventional SCF calculation using approximately 10% additional
   bands;
3. an SPW-SCF calculation.

After the three calculations are completed, the script automatically
reports the numerical accuracy and diagonalization timing information.

Precomputed symmetry files are also provided for the representative
test systems:

```text
S_sg14.mat
S_sg157.mat
S_sg220.mat
```

These files can be used if the MATLAB-Python interface or `pymatgen`
symmetry analysis is unavailable.

### `spw_example/Kpoint`

Contains the finite-$\mathbf{k}$ example used to test the SPW implementation away from the $\Gamma$ point.

The calculation uses the same SPW framework, with symmetry determined from the little group of each irreducible $\mathbf{k}$ point.

Please see the scripts in this directory for the corresponding calculation settings.

### `spw_example/structure_spg`

Contains the crystal-structure setup files used by the example calculations.

For the three representative $\Gamma$-point benchmark systems, the corresponding structure setup files are

```text
test_sg14_natom208_Sb9S4F39_setup.m
test_sg157_natom189_Na21S7Cl_O14F3_2_setup.m
test_sg220_natom156_Cs5Ag4C8IN8_setup.m
```

## Requirements

KSSOLV-SPW requires:

- MATLAB;
- the KSSOLV environment included in this repository;
- Python with `pymatgen` for automatic symmetry analysis.

The Python interface is used to obtain symmetry operations required for

- irreducible Brillouin-zone (IBZ) construction;
- charge-density symmetrization;
- SPW construction.

Users should configure their own Python environment in MATLAB, for example:

```matlab
pyenv(Version="/path/to/python/environment/bin/python");
```

and set the number of MATLAB computational threads according to the available hardware, for example:

```matlab
maxNumCompThreads(48);
```

The value `48` corresponds to the setting used for the calculations in
the associated article and should be changed as appropriate for other
computational environments.

## Getting started

Start MATLAB from the KSSOLV-SPW directory and initialize KSSOLV using

```matlab
KSSOLV_startup
```

The supplied example scripts also call the required KSSOLV initialization routines automatically.

To reproduce one of the representative $\Gamma$-point calculations,
go to

```text
spw_example/Gamma/
```

configure `test_spw.m` for the desired system, and run

```matlab
test_spw
```

Detailed instructions for selecting the structure, symmetry setting,
Python environment, output directory, and alternative precomputed
symmetry data are provided in

```text
spw_example/Reproducing the SPW benchmark examples.md
```

## Reproducing the article calculations

The example calculations use the same principal numerical settings as
those reported in the associated article and Supporting Information,
including:

- a target plane-wave kinetic-energy cutoff of 30 hartree;
- a reduced cutoff of 5 hartree for the pilot calculation;
- a conventional calculation using the target number of bands;
- a conventional calculation using approximately 10% additional bands;
- an SPW calculation using the enlarged band setting.

The detailed computational parameters should be taken from the article
and Supporting Information.

Performance measurements depend on the hardware, MATLAB version,
thread configuration, and Python environment and therefore may differ
from the reported wall-clock times.

## Using KSSOLV-SPW for other systems

The supplied examples are intended primarily to reproduce and validate
the calculations reported in the article.

For other materials, users should provide the appropriate:

- crystal structure;
- space-group information;
- basis and origin transformation required by the crystallographic
  setting;
- $\mathbf{k}$-point sampling;
- plane-wave cutoff;
- number of bands;
- SCF and eigensolver settings.

General KSSOLV input parameters and calculation procedures follow the
standard KSSOLV framework. Please refer to the KSSOLV documentation for
additional information.

## Data availability

The processed data underlying the associated study are provided in the
article and Supporting Information.

Raw calculation log files are available from the corresponding authors
upon reasonable request.

## Citation

If you use KSSOLV-SPW in published work, please cite the associated
KSSOLV-SPW article and the original KSSOLV publication.

The citation information for the KSSOLV-SPW article will be updated
after publication.

## License

KSSOLV-SPW is distributed under the BSD 3-Clause License.

This software is based on and extends KSSOLV. Original KSSOLV copyright
and license notices are retained in the corresponding source files.