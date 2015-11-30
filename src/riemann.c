#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "par.h"
#include "hydro.h"
#include "riemann.h"

int set_riemann(struct parList *pars)
{
    int choice = pars->riemann;
    int err = 0;

    if(choice == 0)
        riemann_flux = &lax_friedrichs_flux;
    else if(choice == 1)
        riemann_flux = &hll_flux;
    else if(choice == 2)
        riemann_flux = &hllc_flux;
    else
    {
        err = 1;
        printf("ERROR - Invalid choice for Riemann Solver: %d\n", choice);
    }

    return err;
}

void lax_friedrichs_flux(double primL[], double primR[], double F[], int nq,
                         double x[2], int dir, struct parList *pars)
{
    double sL, sR, sC, s;
    double U[nq], UL[nq], UR[nq], FL[nq], FR[nq];

    prim2cons(primL, UL, x, 1.0, pars);
    prim2cons(primR, UR, x, 1.0, pars);
    flux(primL, FL, x, dir, pars);
    flux(primR, FR, x, dir, pars);
    wave_speeds(primL, primR, &sL, &sR, &sC, x, dir, pars);

    s = fabs(sL)>fabs(sR) ? fabs(sL) : fabs(sR);

    int q;
    for(q=0; q<nq; q++)
        U[q] = 0.5*(UL[q] + UR[q] - (FR[q] - FL[q])/s);
    for(q=0; q<nq; q++)
        F[q] = 0.5*(FL[q] + FR[q] - s*(UR[q] - UL[q]));
}

void hll_flux(double primL[], double primR[], double F[], int nq,
                double x[2], int dir, struct parList *pars)
{
    double sL, sR, sC;
    double U[nq], UL[nq], UR[nq], FL[nq], FR[nq];

    prim2cons(primL, UL, x, 1.0, pars);
    prim2cons(primR, UR, x, 1.0, pars);
    flux(primL, FL, x, dir, pars);
    flux(primR, FR, x, dir, pars);
    wave_speeds(primL, primR, &sL, &sR, &sC, x, dir, pars);

    int q;
    if(sL > 0)
        for(q=0; q<nq; q++)
        {
            U[q] = UL[q];
            F[q] = FL[q];
        }
    else if(sR < 0)
        for(q=0; q<nq; q++)
        {
            U[q] = UR[q];
            F[q] = FR[q];
        }
    else
        for(q=0; q<nq; q++)
        {
            U[q] = (sR*UR[q] - sL*UL[q] - (FR[q] - FL[q])) / (sR - sL);
            F[q] = (sR*FL[q] - sL*FR[q] + sL*sR*(UR[q] - UL[q])) / (sR - sL);
        }
}

void hllc_flux(double primL[], double primR[], double F[], int nq,
                double x[2], int dir, struct parList *pars)
{
    /*
    double sL, sR, sC;
    double U[nq], UL[nq], UR[nq], FL[nq], FR[nq];

    prim2cons(primL, UL, x, 1.0, pars);
    prim2cons(primR, UR, x, 1.0, pars);
    flux(primL, FL, x, pars);
    flux(primR, FR, x, pars);
    wave_speeds(primL, primR, &sL, &sR, &sC, x, pars);

    int q;
    if(sL > w)
        for(q=0; q<nq; q++)
        {
            U[q] = UL[q];
            F[q] = FL[q];
        }
    else if(sR < w)
        for(q=0; q<nq; q++)
        {
            U[q] = UR[q];
            F[q] = FR[q];
        }
    else if(sL <= w && w <= sC)
    {
        Ustar(primL, U, sL, sC, x, pars);
        for(q=0; q<nq; q++)
            F[q] = FL[q] + sL*(U[q]-UL[q]);
    }
    else if(sC < w && w <= sR)
    {
        Ustar(primR, U, sR, sC, x, pars);
        for(q=0; q<nq; q++)
            F[q] = FR[q] + sR*(U[q]-UR[q]);
    }
    else
    {
        printf("ERROR: wave speeds in HLLC. sL=%.12lg, s=%.12lg, sR=%.12lg\n",
                sL, sC, sR);
    }

    for(q=0; q<nq; q++)
        F[q] -= w*U[q];
    */
}
