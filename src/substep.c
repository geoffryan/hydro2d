#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "par.h"
#include "geom.h"
#include "boundary.h"
#include "hydro.h"
#include "riemann.h"
#include "timestep.h"

//Local Functions

void add_fluxes(struct grid *g, double dt, struct parList *pars);
void add_sources(struct grid *g, double dt, struct parList *pars);

//Definitions

void substep(struct grid *g, double rkfac1, double rkfac2, double dt,
                struct parList *pars)
{
    //Solve Riemann problems.
    add_fluxes(g, dt, pars);

    //Add Sources.
    add_sources(g, dt, pars);

    //Update prims.
    calc_prim(g, pars);

    //Boundary Conditions.
    bc_1L(g, pars);
    bc_1R(g, pars);
    bc_2L(g, pars);
    bc_2R(g, pars);

    //Re-update cons.
    calc_cons(g, pars);
}

//Local Definitions.

void add_fluxes(struct grid *g, double dt, struct parList *pars)
{
    int i, j, q;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int ng11 = g->ng11;
    int ng12 = g->ng12;
    int ng21 = g->ng21;
    int ng22 = g->ng22;
    int d1 = g->d1;
    int d2 = g->d2;
    int nq = g->nq;

    int imin = ng11>0 ? ng11 : 1; 
    int imax = ng12>0 ? nx1-ng12+1 : nx1;
    int jmin = ng21>0 ? ng21 : 1; 
    int jmax = ng22>0 ? nx2-ng22+1 : nx2;

    //X1 Faces.
    for(i = imin; i < imax; i++)
        for(j = ng21; j < nx2-ng22; j++)
        {
            double primL[nq], primR[nq], F[nq];
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i], g->x2[j+1]};
            double x[2];
            geom_CM(xm, xp, x);
            double dA = geom_dA(xm, xp, 0);
            double hn = geom_J(x) / geom_J2(x,0);

            int nL = d1*(i-1)+d2*j;
            int nR = d1*i+d2*j;

            double xL[2] = {g->x1[i-1], g->x2[j]};
            double xR[2] = {g->x1[i+1], g->x2[j+1]};
            double xcL[2], xcR[2];
            geom_CM(xL, xp, xcL);
            geom_CM(xm, xR, xcR);

            for(q=0; q<nq; q++)
            {
                primL[nq] = g->prim[nL+q]
                            + (xp[0]-xcL[0]) * g->prim_grad[2*nL+0+q];
                primR[nq] = g->prim[nR+q]
                            + (xp[0]-xcR[0]) * g->prim_grad[2*nL+0+q];
            }

            reconstruction(g, i, j, 0, primL, primR, pars);

            riemann_flux(primL, primR, F, nq, x, 0, pars);

            for(q=0; q<nq; q++)
            {
                g->cons[d1*(i-1)+d2*j+q] -= F[q] * hn * dA * dt;
                g->cons[d1*  i  +d2*j+q] += F[q] * hn * dA * dt;
            }
        }
    
    //X2 Faces.
    for(i = ng11; i < nx1-ng12; i++)
        for(j = jmin; j < jmax; j++)
        {
            double primL[nq], primR[nq], F[nq];
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i+1], g->x2[j]};
            double x[2];
            geom_CM(xm, xp, x);
            double dA = geom_dA(xm, xp, 1);

            double hn = geom_J(x) / geom_J2(x,1);
            
            reconstruction(g, i, j, 1, primL, primR, pars);
            riemann_flux(primL, primR, F, nq, x, 1, pars);

            for(q=0; q<nq; q++)
            {
                g->cons[d1*i+d2*(j-1)+q] -= F[q] * hn * dA * dt;
                g->cons[d1*i+d2*  j  +q] += F[q] * hn * dA * dt;
            }
        }
}

void add_sources(struct grid *g, double dt, struct parList *pars)
{
    int i, j;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int ng11 = g->ng11;
    int ng12 = g->ng12;
    int ng21 = g->ng21;
    int ng22 = g->ng22;
    int d1 = g->d1;
    int d2 = g->d2;

    for(i=ng11; i<nx1-ng12; i++)
        for(j=ng21; j<nx2-ng22; j++)
        {
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i+1], g->x2[j+1]};

            add_source(&(g->prim[d1*i+d2*j]), &(g->cons[d1*i+d2*j]), 
                        xm, xp, dt, pars);
        }
}

