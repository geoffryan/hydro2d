import math
import sys
import numpy as np
import matplotlib.pyplot as plt
import readCheckpoint as rc

def CM(x1f, x2f, pars):

    geom = pars['Geometry']

    if geom == 1:
        x1 = 0.5*(x1f[1:]+x1f[:-1])
        x2 = 0.5*(x2f[1:]+x2f[:-1])
    elif geom == 2:
        x1 = 0.5*(x1f[1:]+x1f[:-1])
        x2 = 0.5*(x2f[1:]+x2f[:-1])
    elif geom == 3:
        x1p = x1f[1:]
        x1m = x1f[:-1]
        x1 = 2.0*(x1m*x1m+x1m*x1p+x1p*x1p) / (3.0*(x1m+x1p))
        x2 = 0.5*(x2f[1:]+x2f[:-1])
    else:
        x1 = 0.5*(x1f[1:]+x1f[:-1])
        x2 = 0.5*(x2f[1:]+x2f[:-1])

    return x1, x2

def DV(x1f, x2f, pars):

    geom = pars['Geometry']

    if geom == 1:
        dV = 0.5*(x1f[1:]+x1f[:-1])*(x1f[1:]-x1f[:-1])[:,None] \
                * (x2f[1:]-x2f[:-1])[:,None]
    elif geom == 2:
        dV = (np.cos(x1f[:-1])-np.cos(x1f[1:]))[:,None] \
                * (x2f[1:]-x2f[:-1])[None,:]
    elif geom == 3:
        dV = 0.5*(x1f[1:]+x1f[:-1])*(x1f[1:]-x1f[:-1])[:,None] \
                * (x2f[1:]-x2f[:-1])[None,:]
    else:
        dV = (x1f[1:]-x1f[:-1])[:,None] * (x2f[1:]-x2f[:-1])[None,:]

    return dV

def genGrid(x1, x2, pars):
    X1, X2 = np.meshgrid(x1, x2, indexing='ij')
    return X1, X2

def coord2Cart(x1, x2, pars):
    geom = pars['Geometry']

    if geom == 1 or geom == 3:
        x = x1 * np.cos(x2)
        y = x1 * np.sin(x2)
    elif geom == 2:
        x = x1 * np.cos(x2)
        y = x1 * np.sin(x2)
    else:
        x = x1.copy()
        y = x2.copy()

    return x, y

def velCart2Coord(x1, x2, vx, vy, pars):
    geom = pars['Geometry']

    if geom == 1 or geom == 3:
        cosp = np.cos(x2)
        sinp = np.sin(x2)
        vx1 = cosp*vx + sinp*vy
        vx2 = (-sinp*vx + cosp*vy) / x1
    elif geom == 2:
        cosp = np.cos(x2)
        sinp = np.sin(x2)
        vx1 = cosp*vx + sinp*vy
        vx2 = (-sinp*vx + cosp*vy) / x1
    else:
        vx1 = vx.copy()
        vx2 = vy.copy()

    return vx1, vx2

def isentrope1D(X, t, pars):

    gam = pars['GammaLaw']
    rho0 = pars['InitPar1']
    P0 = pars['InitPar2']
    x0 = pars['InitPar3']
    L = pars['InitPar4']
    a = pars['InitPar5']

    K = P0 / math.pow(rho0, gam)
    cs = math.sqrt(gam*P0/rho0)
    J = -2*cs / (gam-1.0)

    rho_min = rho0
    rho_max = rho0*(1.0+a)
    cs_min = math.sqrt(gam * K * math.pow(rho_min, gam-1.0))
    cs_max = math.sqrt(gam * K * math.pow(rho_max, gam-1.0))
    v_min = J + 2.0*cs_min/(gam-1.0)
    v_max = J + 2.0*cs_min/(gam-1.0)

    Xb = X - 0.5*(cs_min+v_min)*t
    Xa = X - 2.0*(cs_max+v_max)*t

    Ntot = np.array(X.shape).prod()
    avgSpace = math.fabs((X.max() - X.min()) / Ntot)
    eps = 1.0e-10
    rho1 = np.ones(X.shape)
    P1 = np.ones(X.shape)
    v1 = np.ones(X.shape)
    cs1 = np.ones(X.shape)

    i = 0
    while np.fabs(Xb-Xa).sum()/Ntot > eps*avgSpace:
        X0 = 0.5*(Xa+Xb)
        dx = (X0-x0)/L
        inInd = np.fabs(dx) < 1.0

        rho1[:] = rho0
        rho1[inInd] += a*rho0*np.power(dx[inInd]*dx[inInd]-1.0, 4)

        P1[:] = K * np.power(rho1, gam)
        cs1[:] = np.sqrt(gam*P1/rho1)
        v1[:] = J + 2.0 * cs1 / (gam-1.0)

        Xt = X0 + (v1+cs1)*t

        overInd = Xt > X
        underInd = Xt < X
        Xb[overInd] = X0[overInd]
        Xa[underInd] = X0[underInd]

        i += 1
        #print i, Xa.min(), Xb.max()

    X0 = 0.5*(Xa+Xb)
    dx = (X0-x0)/L
    inInd = np.fabs(dx) < 1.0

    rho = np.ones(X.shape)
    rho[:] = rho0
    rho[inInd] += a*rho0*np.power(dx[inInd]*dx[inInd]-1.0, 4)

    P = K * np.power(rho, gam)
    cs = np.sqrt(gam*P/rho)
    v = J + 2.0 * cs / (gam-1.0)

    return rho, v, P, X0

def isentrope1D_evolve(X0, t, pars):

    gam = pars['GammaLaw']
    rho0 = pars['InitPar1']
    P0 = pars['InitPar2']
    x0 = pars['InitPar3']
    L = pars['InitPar4']
    a = pars['InitPar5']

    K = P0 / math.pow(rho0, gam)
    cs = math.sqrt(gam*P0/rho0)
    J = -2*cs / (gam-1.0)

    dx = (X0-x0)/L
    inInd = np.fabs(dx) < 1.0

    rho = np.ones(X0.shape)
    rho[:] = rho0
    rho[inInd] += a*rho0*np.power(dx[inInd]*dx[inInd]-1.0, 4)

    P = K * np.power(rho, gam)
    cs = np.sqrt(gam*P/rho)
    v = J + 2.0 * cs / (gam-1.0)

    Xt = X0 + (v+cs)*t

    return rho, v, P, Xt

def isentrope(x1f, x2f, t, pars):

    phi = pars['InitPar6']

    x1, x2 = CM(x1f, x2f, pars)

    X1, X2 = genGrid(x1, x2, pars)
    X,Y = coord2Cart(X1, X2, pars)

    XX = X*math.cos(phi) + Y*math.sin(phi)

    rho, v, P, X0 = isentrope1D(XX, t, pars)
    vx = v*math.cos(phi)
    vy = v*math.sin(phi)
    v1, v2 = velCart2Coord(X, Y, vx, vy, pars)

    return rho, P, v1, v2, XX

def test():
    rho0 = 1.0
    P0 = 1.0
    x0 = 0.2
    L = 0.1
    a = 0.12
    gam = 5.0/3.0
    t = 0.2

    pars = {'InitPar1': rho0, 
            'InitPar2': P0,
            'InitPar3': x0,
            'InitPar4': L,
            'InitPar5': a,
            'GammaLaw': gam}

    X = np.linspace(0.0, 1.0, 1000)

    rho, P, v, X0 = isentrope1D(X, t, pars)
    X00 = np.linspace(X0.min(), X0.max(), 1000)
    rho1, P1, v1, Xt = isentrope1D_evolve(X00, t, pars)

    fig, ax = plt.subplots(2,2)

    ax[0,0].plot(X00, rho1, 'r,')
    ax[0,0].plot(X0, rho, 'k,')
    ax[0,0].plot(X, rho, 'b,')
    ax[0,0].plot(Xt, rho1, 'g,')
    ax[0,0].set_xlabel(r"$\rho$")
    ax[0,1].plot(X00, P1, 'r+')
    ax[0,1].plot(X0, P, 'k+')
    ax[0,1].plot(X, P, 'b+')
    ax[0,1].plot(Xt, P1, 'g+')
    ax[0,1].set_xlabel(r"$P$")
    ax[1,0].plot(X00, v1, 'r+')
    ax[1,0].plot(X0, v, 'k+')
    ax[1,0].plot(X, v, 'b+')
    ax[1,0].plot(Xt, v1, 'g+')
    ax[1,0].set_xlabel(r"$v$")

    plt.show()

def analyze(fname, p=1, plot=True):

    t, x1, x2, prim, cons, pars = rc.readCheckpoint(fname)

    n1 = pars['Nx1']
    n2 = pars['Nx2']
    l1 = pars['X1max'] - pars['X1min']
    l2 = pars['X2max'] - pars['X2min']

    rho_e, P_e, v1_e, v2_e, XX = isentrope(x1, x2, t, pars)

    rho = prim[:,:,0]
    P = prim[:,:,1]
    v1 = prim[:,:,2]
    v2 = prim[:,:,3]

    if plot:
        fig, ax = plt.subplots(2,2)
        ax[0,0].plot(XX, rho_e, 'k+')
        ax[0,0].plot(XX, rho, 'b+')
        ax[0,1].plot(XX, P_e, 'k+')
        ax[0,1].plot(XX, P, 'b+')
        ax[1,0].plot(XX, v1_e, 'k+')
        ax[1,0].plot(XX, v1, 'b+')
        ax[1,1].plot(XX, v1_e, 'k+')
        ax[1,1].plot(XX, v2, 'b+')
    else:
        fig = None

    dV = DV(x1, x2, pars)

    err_rho = math.pow((np.power(np.fabs(rho-rho_e), p)*dV).sum(), 1.0/p)
    err_P = math.pow((np.power(np.fabs(P-P_e), p)*dV).sum(), 1.0/p)
    err_v1 = math.pow((np.power(np.fabs(v1-v1_e), p)*dV).sum(), 1.0/p)
    err_v2 = math.pow((np.power(np.fabs(v2-v2_e), p)*dV).sum(), 1.0/p)

    return n1, n2, l1, l2, err_rho, err_P, err_v1, err_v2, fig

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "Need some checkpoints, bud!"
        test()
        sys.exit()

    
    N1 = []
    N2 = []
    L1 = []
    L2 = []
    ERR_RHO = []
    ERR_P = []
    ERR_V1 = []
    ERR_V2 = []

    for fname in sys.argv[1:]:
        print("Processing {0:s}...".format(fname))
        n1, n2, l1, l2, rho_e, P_e, v1_e, v2_e, fig = analyze(fname)
        N1.append(n1)
        N2.append(n2)
        L1.append(l1)
        L2.append(l2)
        ERR_RHO.append(rho_e)
        ERR_P.append(P_e)
        ERR_V1.append(v1_e)
        ERR_V2.append(v2_e)
        name = ".".join(fname.split("/")[-1].split(".")[:-1])
        fig.savefig(name+".png")
        plt.close(fig)

    print("Error Analysis...")

    N1 = np.array(N1)
    N2 = np.array(N2)
    L1 = np.array(L1)
    L2 = np.array(L2)
    ERR_RHO = np.array(ERR_RHO)
    ERR_P = np.array(ERR_P)
    ERR_V1 = np.array(ERR_V1)
    ERR_V2 = np.array(ERR_V2)

    if len(np.unique(N2)) == 1:
        N = N1.copy()
    elif len(np.unique(N1)) == 1:
        N = N2.copy()
    else:
        N = np.sqrt(N1*N2)

    NN = np.logspace(math.log10(N.min()), math.log10(N.max()), 100)

    fig, ax = plt.subplots(2,2)
    ax[0,0].plot(NN, ERR_RHO[0] * np.power(NN/N[0], -1), ls='--', lw=2.0, 
                    color='grey')
    ax[0,0].plot(NN, ERR_RHO[0] * np.power(NN/N[0], -2), ls='--', lw=2.0, 
                    color='grey')
    ax[0,0].plot(NN, ERR_RHO[0] * np.power(NN/N[0], -3), ls='--', lw=2.0, 
                    color='grey')
    ax[0,0].plot(N, ERR_RHO, 'k+')
    ax[0,1].plot(NN, ERR_P[0] * np.power(NN/N[0], -1), ls='--', lw=2.0, 
                    color='grey')
    ax[0,1].plot(NN, ERR_P[0] * np.power(NN/N[0], -2), ls='--', lw=2.0, 
                    color='grey')
    ax[0,1].plot(NN, ERR_P[0] * np.power(NN/N[0], -3), ls='--', lw=2.0, 
                    color='grey')
    ax[0,1].plot(N, ERR_P, 'k+')
    ax[1,0].plot(NN, ERR_V1[0] * np.power(NN/N[0], -1), ls='--', lw=2.0, 
                    color='grey')
    ax[1,0].plot(NN, ERR_V1[0] * np.power(NN/N[0], -2), ls='--', lw=2.0, 
                    color='grey')
    ax[1,0].plot(NN, ERR_V1[0] * np.power(NN/N[0], -3), ls='--', lw=2.0, 
                    color='grey')
    ax[1,0].plot(N, ERR_V1, 'k+')
    ax[1,1].plot(NN, ERR_V2[0] * np.power(NN/N[0], -1), ls='--', lw=2.0, 
                    color='grey')
    ax[1,1].plot(NN, ERR_V2[0] * np.power(NN/N[0], -2), ls='--', lw=2.0, 
                    color='grey')
    ax[1,1].plot(NN, ERR_V2[0] * np.power(NN/N[0], -3), ls='--', lw=2.0, 
                    color='grey')
    ax[1,1].plot(N, ERR_V2, 'k+')

    ax[0,0].set_xscale('log')
    ax[0,1].set_xscale('log')
    ax[1,0].set_xscale('log')
    ax[1,1].set_xscale('log')
    if (ERR_RHO > 0).all():
        ax[0,0].set_yscale('log')
    if (ERR_P > 0).all():
        ax[0,1].set_yscale('log')
    if (ERR_V1 > 0).all():
        ax[1,0].set_yscale('log')
    if (ERR_V2 > 0).all():
        ax[1,1].set_yscale('log')

    ax[0,0].set_xlabel(r'$\langle N \rangle$')
    ax[0,1].set_xlabel(r'$\langle N \rangle$')
    ax[1,0].set_xlabel(r'$\langle N \rangle$')
    ax[1,1].set_xlabel(r'$\langle N \rangle$')
    ax[0,0].set_ylabel(r'$L_1 (\rho)$')
    ax[0,1].set_ylabel(r'$L_1 (P)$')
    ax[1,0].set_ylabel(r'$L_1 (v^1)$')
    ax[1,1].set_ylabel(r'$L_1 (v^2)$')

    fig.tight_layout()

    fig.savefig("isentrope_convergence.png")




