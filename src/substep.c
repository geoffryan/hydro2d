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
    //Calculate gradients of grid quantities
    calc_grad(g, pars);

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
            double xLm[2] = {g->x1[i-1], g->x2[j]};
            double xLp[2] = {g->x1[i], g->x2[j+1]};
            double xRm[2] = {g->x1[i], g->x2[j]};
            double xRp[2] = {g->x1[i+1], g->x2[j+1]};
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i], g->x2[j+1]};
            double x[2], xL[2], xR[2];
            geom_CM(xLm, xLp, xL);
            geom_CM(xRm, xRp, xR);
            geom_CM2(xm, xp, x, 0);
            double dA = geom_dA(xm, xp, 0);
            double hn = geom_J(x) / geom_J2(x,0);

            int nL = d1*(i-1)+d2*j;
            int nR = d1*i+d2*j;

            for(q=0; q<nq; q++)
            {
                primL[q] = g->prim[nL+q]
                            + (x[0]-xL[0]) * g->prim_grad[2*nL+0*nq+q]
                            + (x[1]-xL[1]) * g->prim_grad[2*nL+1*nq+q];
                primR[q] = g->prim[nR+q]
                            + (x[0]-xR[0]) * g->prim_grad[2*nR+0*nq+q]
                            + (x[1]-xR[1]) * g->prim_grad[2*nR+1*nq+q];
            }

            riemann_flux(primL, primR, F, nq, x, 0, pars);

            for(q=0; q<nq; q++)
            {
                g->cons[nL+q] -= F[q] * hn * dA * dt;
                g->cons[nR+q] += F[q] * hn * dA * dt;
            }

            /*
            printf("face %d: %.12lg %.12lg %.12lg %.12lg\n", i, g->prim[nL+1],
                        primL[1], primR[1], g->prim[nR+1]);
            printf("%d: -%.12lg (%.12lg)\n", i-1, F[2]*hn*dA, F[2]);
            printf("%d: +%.12lg (%.12lg)\n", i, F[2]*hn*dA, F[2]);
            */
        }
    
    //X2 Faces.
    for(i = ng11; i < nx1-ng12; i++)
        for(j = jmin; j < jmax; j++)
        {
            double primL[nq], primR[nq], F[nq];
            double xLm[2] = {g->x1[i], g->x2[j-1]};
            double xLp[2] = {g->x1[i+1], g->x2[j]};
            double xRm[2] = {g->x1[i], g->x2[j]};
            double xRp[2] = {g->x1[i+1], g->x2[j+1]};
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i+1], g->x2[j]};
            double x[2], xL[2], xR[2];
            geom_CM(xLm, xLp, xL);
            geom_CM(xRm, xRp, xR);
            geom_CM2(xm, xp, x, 1);
            
            double dA = geom_dA(xm, xp, 1);
            double hn = geom_J(x) / geom_J2(x,1);
            
            int nL = d1*i+d2*(j-1);
            int nR = d1*i+d2*j;

            for(q=0; q<nq; q++)
            {
                primL[q] = g->prim[nL+q]
                            + (x[0]-xL[0]) * g->prim_grad[2*nL+0*nq+q]
                            + (x[1]-xL[1]) * g->prim_grad[2*nL+1*nq+q];
                primR[q] = g->prim[nR+q]
                            + (x[0]-xR[0]) * g->prim_grad[2*nR+0*nq+q]
                            + (x[1]-xR[1]) * g->prim_grad[2*nR+1*nq+q];
            }

            riemann_flux(primL, primR, F, nq, x, 1, pars);

            for(q=0; q<nq; q++)
            {
                g->cons[nL+q] -= F[q] * hn * dA * dt;
                g->cons[nR+q] += F[q] * hn * dA * dt;
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
    int nq = g->nq;

    for(i=ng11; i<nx1-ng12; i++)
        for(j=ng21; j<nx2-ng22; j++)
        {
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i+1], g->x2[j+1]};

            //printf("%d:", i);

            add_source(&(g->prim[d1*i+d2*j]), &(g->cons[d1*i+d2*j]), 
                        &(g->prim_grad[2*(d1*i+d2*j)   ]), 
                        &(g->prim_grad[2*(d1*i+d2*j)+nq]), xm, xp, dt, pars);
        }
}

