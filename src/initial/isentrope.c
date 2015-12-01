#include <stdlib.h>
#include <math.h>
#include "../geom.h"
#include "../initial.h"
#include "../par.h"
#include "../hydro/newt2.h"

// Initial condition for an isentropic wave.  See Section 4.6
// of the RAM paper (Zhang & MacFadyen 2006) for details.

void initial_isentrope(double *prim, double x[2], int nq, struct parList *pars)
{
    double rho0 = pars->initPar1; // Reference density
    double P0 = pars->initPar2;   // Reference pressure
    double x0 = pars->initPar3;   // Pulse location
    double L = pars->initPar4;    // Pulse width
    double a = pars->initPar5;    // Pulse strength
    double phi = pars->initPar6;    // Pulse Direction

    double rho, P, v, dx, K, GAM, cs, J;

    GAM = pars->gammalaw;
    K = P0 / pow(rho0, GAM);
    dx = (x[0]*cos(phi)+x[1]*sin(phi)-x0)/L;
    cs = sqrt(GAM*P0/rho0);

    J = -2.0 * cs / (GAM-1.0); 

    rho = rho0;
    P = P0;
    v = 0.0;

    if(fabs(dx) < 1.0)
    {
        rho += a * rho0 * pow(dx*dx-1.0, 4);
        P = K * pow(rho, GAM);
        cs = sqrt(GAM*P/rho); 
        v = J +  2.0 * cs / (GAM-1.0);
    }

    double vx1, vx2;

    vx1 = v*cos(phi);
    vx2 = v*sin(phi);

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[VX1] = vx1;
    prim[VX2] = vx2;

    int q;
    for(q=4; q<nq; q++)
        prim[q] = dx>0 ? 0.0 : 1.0;
}
