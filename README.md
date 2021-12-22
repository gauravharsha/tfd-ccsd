# tfd-ccsd
### Thermofield dynmaics based coupled cluster with singles and doubles

Python package for computing grand-canonical ensemble properties of many-body fermionic systems at finite-temperatures using the CCSD approximation. The theory and equations are discussed in the following paper

> Thermofield Theory for Finite-Temperature Coupled Cluster, G. Harsha, T. M. Henderson, and G. E. Scuseria, [J. Chem. Theory Comput. 15, 6127-6136 (2019)](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00744) [arXiv](https://arxiv.org/abs/1907.11286)

This code performs the imaginary time evolution in the interaction picture, i.e., using the covariant coupled cluster ansatz. A fixed-reference version can be formulated by tweaking the evolution equations here and there but is not explicitly included in this repository.

## Requirements

The package requires Python 3.7 or higher. Computationally intensive modules and functions that drives the imaginary time and chemical potential evolution are written in Fortran. The equations were generated using a modified [drudge](https://github.com/tschijnmo/drudge) class called `ThermofieldDrudge`. Drudge scripts are not included in this repository, but the interested user may check out the [thermal-ci](https://github.com/gauravharsha/thermal-ci) repository for an example.

## Installation

If you already have the packages listed in `requirements.txt`, you can simply install this package by
```bashscript
python setup.py install
```

## Usage

1. Prepare the one- and two-electron integrals in the basis in which mean-field Hamiltonian is diagonal. We will also need these mean-field eigenvalues. In the current implementation, the thermal-ci package expects these integrals in the spin-orbital basis, stored as an HDF (.h5) file with the keys `h1`, `eri` and `eigs` respectively.

   `h1`   :   2D array of size (N,N)
   
   `eri`  :   4D array of size (N,N,N,N)
   
   `eigs` :   1D array of size (N)
   
   The integrals can be generated using standard Quantum Chemistry packages such as [pyscf](https://github.com/pyscf/pyscf). For example, for 2-site Hubbard model with U/t = 1, we can have a data file named `hub_2s_u1_data.h5`.

   Each of these objects in the HDF5 integrals file should have the same list of attributes.

2. With the integrals' file ready, the details about the system and the thermal evolution needs to be specified in the `Input` file. As an example, see the the Input file in the `scripts` subdirectory.

3. Scripts to perform the imaginary time evolution and compute internal energies within this thermal CCSD approximation are included in the  `scripts` subdirectory.

   `main_fixedN.py`     :   Python script for imaginary time evolution using covariant CCSD while maintaining a fixed number of particles.

   `main_fixedMu.py`    :   Python script for imaginary time evolution using covariant CCSD for a set value of the chemical potential.
