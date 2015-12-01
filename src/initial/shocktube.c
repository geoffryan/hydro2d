#include <stdlib.h>
#include <math.h>
#include "../geom.h"
#include "../initial.h"
#include "../par.h"
#include "../hydro/newt2.h"

// Initial condition for a generic shocktube

void initial_shocktube(double *prim, double x[2], int nq, struct parList *pars)
{
    double rhoL = pars->initPar1;   // Left density
    double PL = pars->initPar2;     // Left pressure
    double vL = pars->initPar3;     // Left velocity
    double rhoR = pars->initPar4;   // Right density
    double PR = pars->initPar5;     // Right pressure
    double vR = pars->initPar6;     // Right velocity
    double x0 = pars->initPar7;     // Shock position
    double phi = pars->initPar8;     // Shock direction

    double rho, P, v, dx;

    dx = x[0]*cos(phi)+x[1]*sin(phi)-x0;

    if(dx > 0)
    {
        rho = rhoR;
        P = PR;
        v = vR;
    }
    else
    {
        rho = rhoL;
        P = PL;
        v = vL;
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
