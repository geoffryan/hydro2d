import sys
import numpy as np
import matplotlib.pyplot as plt
import readCheckpoint as rc

def plot2DSingle(fig, ax, X1, X2, dat, title=None):
    im = ax.pcolormesh(X1, X2, dat, cmap=plt.cm.viridis)
    ax.set_xlabel(r'$x^1$')
    ax.set_ylabel(r'$x^2$')
    fig.colorbar(im, ax=ax)
    if title is not None:
        ax.set_title(title)


def plotAll2D(filename):

    dat = rc.readCheckpoint(filename)
    print("{0:s}: {1:s}".format(filename, rc.summary(*dat)))

    t = dat[0]
    x1 = dat[1]
    x2 = dat[2]
    prim = dat[3]
    pars = dat[5]

    X1, X2 = np.meshgrid(x1, x2, indexing='ij')

    fig, ax = plt.subplots(2,2,figsize=(12,9))

    plot2DSingle(fig, ax[0,0], X1, X2, prim[:,:,0])
    plot2DSingle(fig, ax[0,1], X1, X2, prim[:,:,1])
    plot2DSingle(fig, ax[1,0], X1, X2, prim[:,:,2])
    plot2DSingle(fig, ax[1,1], X1, X2, prim[:,:,3])

    fig.suptitle("t = {0:.3f} ({1:s})".format(t, pars["GitHash"]))

    root = ".".join(filename.split("/")[-1].split(".")[:-1])
    plotname = "{0:s}_{1:s}.png".format("plot2d", root)

    print('    Saving: {0:s}'.format(plotname))
    fig.savefig(plotname)

    plt.close(fig)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: python plotBasic.py <checkpoint.h5 ...>")
        print("    Creates 2D plots of primitive quantities.")
        sys.exit()

    for file in sys.argv[1:]:
        plotAll2D(file)


