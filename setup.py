from numpy.distutils.core import setup, Extension


thermalcisd = Extension(
    'tfdccsd.ThermalCISD',
    sources=['tfdccsd/fort_src/ThermalCISD.f90', ],
    libraries=['lapack', 'blas'],
    extra_f90_compile_args=["-O3", "-fopenmp"],
    extra_link_args=["-lgomp"]
)

thermalccsd = Extension(
    'tfdccsd.ThermalCCSD',
    sources=['tfdccsd/fort_src/ThermalCCSD.f90', ],
    libraries=['lapack', 'blas'],
    extra_f90_compile_args=["-O3", "-fopenmp"],
    extra_link_args=["-lgomp"]
)

expvals = Extension(
    'tfdccsd.ExpVals',
    sources=['tfdccsd/fort_src/ExpVals.f90', ],
    libraries=['lapack', 'blas'],
    extra_f90_compile_args=["-O3", "-fopenmp"],
    extra_link_args=["-lgomp"]
)

setup(
   name='tfdccsd',
   version='1.2',
   description='A package for calculating finite temperature coupled cluster\
       properties of many-body fermionic systems in the\
       grand canonical ensemble',
   packages=['tfdccsd'],
   ext_modules=[thermalccsd, thermalcisd, expvals]
)
