#!/bin/csh

# Switch to fort_src directory
cd fort_src

# Use f2py to form the python modules for the Covariant CISD
# f2py -c --verbose --opt='-O4' -m ThermalCISD ThermalCISD.f90 --f90flags="-fopenmp" -lgomp
f2py -c --verbose --opt='-O4' -m ThermalCISD ThermalCISD.f90 --fcompiler=pg --f90flags="-openmp" -lgomp -lblas

# Use f2py to form the python modules for the Covariant CCSD
# f2py -c --verbose --opt='-O3' -m ThermalCCSD ThermalCCSD.f90 --f90flags="-fopenmp" -lgomp -lblas
f2py -c --verbose --opt='-O4' -m ThermalCCSD ThermalCCSD.f90 --fcompiler=pg --f90flags="-openmp" -lgomp -lblas

# Use f2py to form the python modules for the Expectation Values
# f2py -c --verbose --opt='-O4' -m ExpVals ExpVals.f90 --f90flags="-fopenmp" -lgomp
f2py -c --verbose --opt='-O4' -m ExpVals ExpVals.f90 --fcompiler=pg --f90flags="-openmp" -lgomp

# Get out to main directory
cd ..
