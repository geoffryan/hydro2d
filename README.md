# 2D Hydrodynamics

This is a test code to play around with geometry implementations in
hydrodynamic codes.  It solves the two dimensional adiabatic Euler equations
in arbitrary geometety, in theory.

## Installation:

Copy `Makefile.in.template` to `Makefile.in`.  Edit `Makefile.in` to point to
your local HDF5 installation and specify your C compiler.  Then just run 
`make`. The executable is installed in `bin/`.

`Makefile.in` is ignored by the repo, it is machine specific, so feel free to 
edit it at any time.

To automatically copy the executable, parfiles, and visualization scripts into
another directory set the `INSTALL_DIR` directory in `Makefile.in` and run 
`make install`.

## Running:

    $ bin/hydro2d myparfile.par

Sample parfiles are included in `parfiles/`.

## Visualization & Utilities

Some python scripts for plotting and examining output are included in `vis/`.
The `h5py` module is required for reading the HDF5 format output, and `numpy`
and `matplotlib` are used for plotting and data manipulation.


