# Reproducing the SPW benchmark examples

The script `test_spw.m` is provided to reproduce the representative SPW calculations reported in the associated article.

Three example crystalline systems are included. For each selected system, `test_spw` performs three calculations sequentially:

1. a conventional SCF calculation using the target number of bands;
2. a conventional SCF calculation using approximately 10% additional bands;
3. an SPW-SCF calculation using the enlarged band setting.

After all three calculations are completed, the script automatically reports the eigenvalue and charge-density differences and the diagonalization timing information.

## Requirements

`test_spw.m` requires a working KSSOLV-SPW installation together with the dependencies required by KSSOLV.

A Python environment containing `pymatgen` is also required for the automatic symmetry analysis. The Python interface is used to identify the crystal symmetry operations required by the irreducible Brillouin-zone (IBZ) construction, charge-density symmetrization, and SPW construction.

The calculations reported in the article used:

- target plane-wave kinetic-energy cutoff: 30 hartree;
- reduced-cutoff pilot calculation: 5 hartree.

Other numerical settings are consistent with those described in the article and Supporting Information.

For general KSSOLV input parameters and usage, please refer to the KSSOLV documentation.

## 1. Configure the computational environment

Before running the calculation, modify the MATLAB thread limit and Python environment according to your local computational environment.

For example:

```matlab
maxNumCompThreads(48);

pyenv(Version="/path/to/python/environment/bin/python");
```

The selected Python environment must contain `pymatgen`.

The number of MATLAB computational threads does not need to be 48. Users should select an appropriate value for their own hardware.

## 2. Select the output directory

Set `folder_path` to the directory in which the calculation log and related output files should be stored.

For example:

```
folder_path = fullfile(pwd,'sg157-C3v');
```

Use a different output directory for each test system to avoid overwriting previous results.

## 3. Select one test material

Three representative systems from the article are provided.

Uncomment exactly one of the following structure inputs and comment out the others:

```
% Space group 14
test_sg14_natom208_Sb9S4F39_setup;

% Space group 157
test_sg157_natom189_Na21S7Cl_O14F3_2_setup;

% Space group 220
test_sg220_natom156_Cs5Ag4C8IN8_setup;
```

Only one structure setup should be active for each run.

## 4. Configure the corresponding space-group information

The crystal symmetry is first obtained automatically by

```
S = Sym(cry);
```

The corresponding Bilbao space-group information must then be added using `Sbilbaoinfo`. The space-group number and basis/origin transformation must match the selected structure.

Use the following settings for the three article examples:

| Example          | Space group | Structure input                              | Symmetry setting                     |
| ---------------- | ----------- | -------------------------------------------- | ------------------------------------ |
| Sb9S4F39         | 14          | `test_sg14_natom208_Sb9S4F39_setup`          | `S = Sbilbaoinfo(S,14,eye(3),0*p);`  |
| Na21S7Cl(O14F3)2 | 157         | `test_sg157_natom189_Na21S7Cl_O14F3_2_setup` | `S = Sbilbaoinfo(S,157,eye(3),0*p);` |
| Cs5Ag4C8IN8      | 220         | `test_sg220_natom156_Cs5Ag4C8IN8_setup`      | `S = Sbilbaoinfo(S,220,P_bcc,0*p);`  |

The matrices such as `P_bcc` specify the required basis transformation for the corresponding crystallographic setting. The last argument specifies the origin transformation.

For example, to reproduce the space-group-157 system:

```
test_sg157_natom189_Na21S7Cl_O14F3_2_setup;

...

S = Sym(cry);
S = Sbilbaoinfo(S,157,eye(3),0*p);
cry = genIBZ(cry,S);
```

Make sure that the selected structure input and the `Sbilbaoinfo` setting always correspond to the same material.

The corresponding settings in `test_spw.m` are:

```
% S = Sbilbaoinfo(S,14,eye(3),0*p);
% S = Sbilbaoinfo(S,157,eye(3),0*p);
% S = Sbilbaoinfo(S,220,P_bcc,0*p);
```

## 5. Alternative symmetry input

On some systems, the MATLAB-Python interface or `pymatgen` may not initialize correctly.

Precomputed symmetry information is provided as an alternative for the example calculations.

In this case, replace

```
S = Sym(cry);
```

with the corresponding supplied symmetry file.

For example:

```
load('S_sg157.mat','S');
```

Then continue with the appropriate `Sbilbaoinfo` and `genIBZ` calls:

```
S = Sbilbaoinfo(S,157,eye(3),0*p);
cry = genIBZ(cry,S);
```

This fallback allows the example calculation to be reproduced without performing the automatic `pymatgen` symmetry search.

## 6. Run the calculation

After completing the settings above, run

```
test_spw
```

from MATLAB.

The script performs the three SCF calculations automatically in the following order:

```
Part 1: Conventional SCF
        Target number of bands

Part 2: Conventional SCF with extra bands
        Approximately 10% additional bands

Part 3: SPW-SCF
        SPW calculation with the enlarged band setting
```

No manual restart is required between these three calculations.

The production SCF convergence tolerance, reduced-cutoff tolerance, eigensolver settings, band enlargement, and other parameters are set inside `test_spw.m` consistently with the calculations reported in the article.

## 7. Output and accuracy report

After all three calculations are finished, `test_spw.m` automatically compares the calculations and prints:

- Kohn-Sham eigenvalue errors;
- charge-density errors;
- diagonalization timing/profile information for the conventional target-band calculation;
- diagonalization timing/profile information for the conventional extra-band calculation;
- diagonalization timing/profile information for the SPW calculation.

The main text output is written to a timestamped log file in the selected output directory, with a name of the form

```
matlab_accuracy_<eigensolver>_<timestamp>.log
```

For example:

```
matlab_accuracy_davidson_qe_20260808_120000.log
```

The density data used for the accuracy comparison are stored temporarily during the calculation and removed after the final report has been generated.

MATLAB profiling can additionally be enabled by setting

```
useprofile = 1;
```

which generates profiling data for the individual calculation steps.

## Notes

The example settings are intended to reproduce the calculations reported in the associated article.

Performance will depend on the hardware, MATLAB configuration, number of computational threads, and software environment.

For calculations on other materials, users should modify the crystal input, space-group information, basis/origin transformation, k-point sampling, plane-wave cutoff, and other KSSOLV parameters as needed.

Please refer to the KSSOLV documentation for general input preparation and calculation settings.
\```