#include <stdlib.h>
#include "../geom.h"
#include "../initial.h"
#include "../par.h"
#include "../hydro/newt2.h"

void initial_uniform(double *prim, double x[2], int nq, struct parList *par)
{
    double rho0 = par->initPar1;
    double P0 = par->initPar2;
    double vx1 = par->initPar3;
    double vx2 = par->initPar4;

    prim[RHO] = rho0;
    prim[PPP] = P0;
    prim[VX1] = vx1;
    prim[VX2] = vx2;

    int q;
    for(q=NC; q<nq; q++)
        prim[q] = 0.0;
}
