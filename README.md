# KSSOLV-SPW

KSSOLV-SPW implements the symmetrized plane-wave (SPW) framework in
KSSOLV to accelerate iterative diagonalization in plane-wave density
functional theory using space-group symmetry.

## Overview

The SPW method decomposes the plane-wave eigenproblem into independent
symmetry-adapted subproblems associated with the small representations
of the relevant little group. The implementation retains the underlying
KSSOLV block-Davidson eigensolver and conventional Hamiltonian
application.

## Requirements

- MATLAB R2022a or later
- KSSOLV
- Pymatgen
- Additional dependencies listed in the source code

## Usage

Instructions and example calculations will be provided here.

## Reproducibility data

The benchmark data reported in the associated publication are described
in the article and Supporting Information. Additional raw calculation
log files are available from the corresponding authors upon reasonable
request.

## Citation

If you use KSSOLV-SPW, please cite:

[To be submitted]

and the original KSSOLV publication:

[KSSOLV citation]
Jiao, S.; Zhang, Z.; Wu, K.; Wan, L.; Ma, H.; Li, J.; Chen, S.; Qin, X.; Liu, J.;
Ding, Z.; Yang, J.; Li, Y.; Hu, W.; Lin, L.; Yang, C. KSSOLV 2.0: An efficient MATLAB
toolbox for solving the Kohn-Sham equations with plane-wave basis set. Comput. Phys.
Commun. 2022, 279, 108424.

## License

KSSOLV-SPW is distributed under the BSD 3-Clause License.

This software is based on and extends KSSOLV Original KSSOLV
copyright and license notices are retained in the corresponding source
files.

## Contact

For questions about KSSOLV-SPW, please contact:

WeiHu whuustc@ustc.edu.cn
