#include <stdlib.h>
#include "../boundary.h"
#include "../geom.h"
#include "../grid.h"
#include "../initial.h"
#include "../par.h"

// Fixed (dirichlet) boundary condition.

void bc_1L_fixed(struct grid *g, struct parList *par)
{
    int nq = g->nq;
    int nx2 = g->nx2;
    int ng11 = g->ng11;

    int i, j;

    for(i=0; i<ng11; i++)
        for(j=0; j<nx2; j++)
        {
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i+1], g->x2[j+1]};
            double x[2];
            geom_CM(xm, xp, x);
            initial_value(&(g->prim[nq*(nx2*i+j)]), x, nq, par);
        }
}

void bc_1R_fixed(struct grid *g, struct parList *par)
{
    int nq = g->nq;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int ng12 = g->ng12;

    int i, j;

    for(i=nx1-ng12; i<nx1; i++)
        for(j=0; j<nx2; j++)
        {
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i+1], g->x2[j+1]};
            double x[2];
            geom_CM(xm, xp, x);
            initial_value(&(g->prim[nq*(nx2*i+j)]), x, nq, par);
        }
}

void bc_2L_fixed(struct grid *g, struct parList *par)
{
    int nq = g->nq;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int ng21 = g->ng21;

    int i, j;

    for(i=0; i<nx1; i++)
        for(j=0; j<ng21; j++)
        {
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i+1], g->x2[j+1]};
            double x[2];
            geom_CM(xm, xp, x);
            initial_value(&(g->prim[nq*(nx2*i+j)]), x, nq, par);
        }
}

void bc_2R_fixed(struct grid *g, struct parList *par)
{
    int nq = g->nq;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int ng22 = g->ng22;

    int i, j;

    for(i=0; i<nx1; i++)
        for(j=nx2-ng22; j<nx2; j++)
        {
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i+1], g->x2[j+1]};
            double x[2];
            geom_CM(xm, xp, x);
            initial_value(&(g->prim[nq*(nx2*i+j)]), x, nq, par);
        }
}
