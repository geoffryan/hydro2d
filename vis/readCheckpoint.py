import sys
import numpy as np
import h5py as h5

# Reads checkpoints created by hydro2d.

def readCheckpoint(filename):

    f = h5.File(filename, "r")
    dat = f['data']

    X1 = dat['x1'][...]
    X2 = dat['x2'][...]
    prim = dat['prim'][...]
    cons = dat['cons'][...]
    t = dat['t'][0]

    pars = {}
    for key in dat.attrs:
        val = dat.attrs[key]
        try:
            if len(val.shape) == 1 and val.shape[0] == 1:
                pars[key] = val[0]
            else:
                pars[key] = val
        except:
            pars[key] = val

    f.close()

    return t, X1, X2, prim, cons, pars

def summary(t, X1, X2, prim, cons, pars):
    #Return simple summary of checkpoint data.
    sum = "t={0:f} ({1:d}x{2:d}) nq={3:d} {4:s}".format(
            t, X1.shape[0]-1, X2.shape[0]-1, prim.shape[2], 
            pars["GitHash"])
    return sum

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: python readCheckpoint.py <checkpoint.h5 ...>")
        sys.exit()

    for file in sys.argv[1:]:
        data = readCheckpoint(file)
        print("{0:s}: {1:s}".format(file, summary(*data)))

