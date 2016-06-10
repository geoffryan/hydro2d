#include <stdlib.h>
#include <math.h>
#include "../geom.h"
#include "../initial.h"
#include "../par.h"
#include "../hydro/newt2.h"

void initial_bomb(double *prim, double x[2], int nq, struct parList *par)
{
    double rho0 = par->initPar1;
    double P0 = par->initPar2;
    double rho1 = par->initPar3;
    double P1 = par->initPar4;
    double x1a = par->initPar5;
    double x1b = par->initPar6;
    double x2a = par->initPar7;
    double x2b = par->initPar8;

    double rho, P;
    if(x[0] > x1a && x[0] < x1b && x[1] > x2a && x[1] < x2b)
    {
        rho = rho1;
        P = P1;
    }
    else
    {
        rho = rho0;
        P = P0;
    }

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[VX1] = 0.0;
    prim[VX2] = 0.0;

    int q;
    for(q=NC; q<nq; q++)
        prim[q] = 0.0;
}
