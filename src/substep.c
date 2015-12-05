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
    int i, j, q;
    /*
    printf("Beginning.\n");
    for(i=0; i<g->nx1; i++)
        for(j=0; j<g->nx2; j++)
            for(q=0; q<g->nq; q++)
                if(g->prim[(g->nq)*((g->nx2)*i+j)+q] 
                        != g->prim[(g->nq)*((g->nx2)*i+j)+q]
                    || g->cons[(g->nq)*((g->nx2)*i+j)+q] 
                        != g->cons[(g->nq)*((g->nx2)*i+j)+q])
                    printf("Cell %d %d val %d is NaN\n", i, j, q);
    */

    //Solve Riemann problems.
    add_fluxes(g, dt, pars);

    /*
    printf("Fluxes added.\n");
    for(i=0; i<g->nx1; i++)
        for(j=0; j<g->nx2; j++)
            for(q=0; q<g->nq; q++)
                if(g->cons[(g->nq)*((g->nx2)*i+j)+q] 
                        != g->cons[(g->nq)*((g->nx2)*i+j)+q])
                    printf("Cell %d %d val %d is NaN\n", i, j, q);
    */

    //Add Sources.
    add_sources(g, dt, pars);

    /*
    printf("Sources added.\n");
    for(i=0; i<g->nx1; i++)
        for(j=0; j<g->nx2; j++)
            for(q=0; q<g->nq; q++)
                if(g->cons[(g->nq)*((g->nx2)*i+j)+q] 
                        != g->cons[(g->nq)*((g->nx2)*i+j)+q])
                    printf("Cell %d %d val %d is NaN\n", i, j, q);
    */

    //Update prims.
    calc_prim(g, pars);
    
    /*
    printf("prims recalculated.\n");
    for(i=0; i<g->nx1; i++)
        for(j=0; j<g->nx2; j++)
            for(q=0; q<g->nq; q++)
                if(g->prim[(g->nq)*((g->nx2)*i+j)+q] 
                        != g->prim[(g->nq)*((g->nx2)*i+j)+q])
                    printf("Cell %d %d val %d is NaN\n", i, j, q);
    */

    //Boundary Conditions.
    bc_1L(g, pars);
    bc_1R(g, pars);
    bc_2L(g, pars);
    bc_2R(g, pars);
    
    /*
    printf("Boundaries applied.\n");
    for(i=0; i<g->nx1; i++)
        for(j=0; j<g->nx2; j++)
            for(q=0; q<g->nq; q++)
                if(g->prim[(g->nq)*((g->nx2)*i+j)+q] 
                        != g->prim[(g->nq)*((g->nx2)*i+j)+q])
                    printf("Cell %d %d val %d is NaN\n", i, j, q);
    */

    //Re-update cons.
    calc_cons(g, pars);
    
    /*
    printf("conss recalculated.\n");
    for(i=0; i<g->nx1; i++)
        for(j=0; j<g->nx2; j++)
            for(q=0; q<g->nq; q++)
                if(g->cons[(g->nq)*((g->nx2)*i+j)+q] 
                        != g->cons[(g->nq)*((g->nx2)*i+j)+q])
                    printf("Cell %d %d val %d is NaN\n", i, j, q);
    */
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

            reconstruction(g, i, j, 0, primL, primR, pars);
            /*
            for(q=0; q<nq; q++)
                if(primL[q] != primL[q])
                    printf("x1 face %d %d PL[%d] is NaN\n", i, j, q);

            for(q=0; q<nq; q++)
                if(primR[q] != primR[q])
                    printf("x1 face %d %d PR[%d] is NaN\n", i, j, q);

            if(i==1 && j ==2)
            {
                printf("primL:");
                for(q=0; q<nq; q++)
                    printf(" %lg", primL[q]);
                printf("\nprimR:");
                for(q=0; q<nq; q++)
                    printf(" %lg", primR[q]);
                printf("\n");
            }
            */

            riemann_flux(primL, primR, F, nq, x, 0, pars);

            /*
            for(q=0; q<nq; q++)
                if(F[q] != F[q])
                    printf("x1 face %d %d F[%d] is NaN\n", i, j, q);
            */

            for(q=0; q<nq; q++)
            {
                g->cons[nq*(nx2*(i-1)+j)+q] -= F[q] * hn * dA * dt;
                g->cons[nq*(nx2* i   +j)+q] += F[q] * hn * dA * dt;
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

            /*
            for(q=0; q<nq; q++)
                if(F[q] != F[q])
                    printf("x2 face %d %d F[%d] is NaN\n", i, j, q);
            */

            for(q=0; q<nq; q++)
            {
                g->cons[nq*(nx2*i+j-1)+q] -= F[q] * hn * dA * dt;
                g->cons[nq*(nx2*i+ j )+q] += F[q] * hn * dA * dt;
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
    int nq = g->nq;

    for(i=ng11; i<nx1-ng12; i++)
        for(j=ng21; j<nx2-ng22; j++)
        {
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i+1], g->x2[j+1]};
            double x[2];
            geom_CM(xm,xp,x);
            double dV = geom_dV(xm,xp);

            add_source(&(g->prim[nq*(nx2*i+j)]), &(g->cons[nq*(nx2*i+j)]), 
                        x, dV*dt, pars);
        }
}

