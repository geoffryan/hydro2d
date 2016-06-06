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
    cs = sqrt(GAM*P0/rho0);

    J = -2.0 * cs / (GAM-1.0); 

    if(pars->geom == 1 || pars->geom == 3)
        dx = (x[0]*cos(x[1])*cos(phi)+x[0]*sin(x[1])*sin(phi)-x0)/L;
    else
        dx = (x[0]*cos(phi)+x[1]*sin(phi)-x0)/L;
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

    double vx1, vx2, vx, vy;
    
    vx = v*cos(phi);
    vy = v*sin(phi);

    if(pars->geom == 1 || pars->geom == 3)
    {
        vx1 = cos(x[1])*vx + sin(x[1])*vy;
        vx2 = (-sin(x[1])*vx + cos(x[1])*vy) / x[0];
    }
    else
    {
        vx1 = vx;
        vx2 = vy;
    }

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[VX1] = vx1;
    prim[VX2] = vx2;

    int q;
    for(q=4; q<nq; q++)
        prim[q] = dx>0 ? 0.0 : 1.0;
}
