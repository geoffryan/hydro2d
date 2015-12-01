#include <stdlib.h>
#include <math.h>
#include "../geom.h"
#include "../initial.h"
#include "../par.h"
#include "../hydro/newt2.h"

void initial_uniform(double *prim, double x[2], int nq, struct parList *par)
{
    double rho0 = par->initPar1;
    double P0 = par->initPar2;
    double vx = par->initPar3;
    double vy = par->initPar4;
    double vx1, vx2;

    if(par->geom == 1)
    {
        double r = x[0];
        double phi = x[1];
        vx1 = cos(phi)*vx + sin(phi)*vy;
        vx2 = (-sin(phi)*vx + cos(phi)*vy) / r;
    }
    else
    {
        vx1 = vx;
        vx2 = vy;
    }

    prim[RHO] = rho0;
    prim[PPP] = P0;
    prim[VX1] = vx1;
    prim[VX2] = vx2;

    int q;
    for(q=NC; q<nq; q++)
        prim[q] = 0.0;
}
