# 2D Hydrodynamics

This is a test code to play around with geometry implementations in
hydrodynamic codes.

## Installation:

Copy `Makefile.in.template` to `Makefile.in`.  Edit `Makefile.in` to point to
your local HDF5 installation and specify your C compiler.  Then just run 
`make`. The executable is installed in `bin/`.

`Makefile.in` is ignored by the repo, it is machine specific, so feel free to 
edit it at any time.

## Running:

    $ bin/hydro2d myparfile.par

Sample parfiles are included in `parfiles/`.


