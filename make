#!/bin/csh

# Switch to fort_src directory
cd fort_src

# Use f2py to form the python modules for the Covariant CISD
f2py -c --verbose --opt='-O3' -m ThermalCISD ThermalCISD.f90 --f90flags="-fopenmp" -lgomp

# Use f2py to form the python modules for the Covariant CCSD
f2py -c --verbose --opt='-O3' -m ThermalCCSD ThermalCCSD.f90 --f90flags="-fopenmp" -lgomp

# Use f2py to form the python modules for the Expectation Values
f2py -c --verbose --opt='-O3' -m ExpVals ExpVals.f90 --f90flags="-fopenmp" -lgomp

# Get out to main directory
cd ..
