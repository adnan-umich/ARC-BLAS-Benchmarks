# HLIBcore

*HLIBcore* contains basic type definitions and function implementations from
[HLIBpro](https://hlibpro.com) needed by [libHLR](http://libhlr.org), e.g., matrix types, clustering
algorithms and input/output functions. 

The API is a fully compatible subset of the API provided by *HLIBpro*, i.e., you can directly
replace *HLIBcore* with *HLIBpro*. However, all the arithmetic functions for H-matrices and even
basic matrix construction is missing.  For this, please use the functions provided by *libHLR*.

## Prerequisites

HLIBcore is tested on Linux and MacOS systems with GCC (>= v5), Clang (>= v4) and Intel
C++ compiler (>= v17). As a general rule: C++14 is required. Also a Fortran compiler is
needed, e.g., `gfortran` or `f77`.

Furthermore, the [Boost](https://boost.org) C++ library and BLAS/LAPACK are needed for HLIBcore.

For compilation, the [SCons](https://scons.org) construction tool is used, while the
configuration utility is based on Python 3.

Optional libraries are

- [Intel MKL](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html)
- [Mongoose](https://github.com/scottkolo/mongoose)
- [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
- [Scotch](https://www.labri.fr/perso/pelegrin/scotch/)
- [GNU Scientific Library](https://www.gnu.org/software/gsl/)
- [HDF5](https://www.hdfgroup.org)

These libraries can normally be installed via the standard packet management systems of
your operating system.

### Debian/Ubuntu

~~~
> apt install gcc g++ gfortran scons python3
> apt install libboost-all-dev liblapack-dev libsuitesparse-dev libmetis-dev libscotch-dev libgsl-dev libhdf5-dev
~~~

## Installation

Download an archived release or clone the repository:

~~~
> git clone https://gitlab.mis.mpg.de/rok/hlibcore
~~~

Afterwards, call the configuration tool

~~~
> ./configure
~~~

If configuration went without problems, you can start compilation

~~~
> scons
~~~
