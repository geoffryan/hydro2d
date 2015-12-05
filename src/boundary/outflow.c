#include <stdlib.h>
#include "../boundary.h"
#include "../grid.h"
#include "../par.h"

// Outflow (Zero coordinate gradient) boundary.

void bc_1L_outflow(struct grid *g, struct parList *par)
{
    int nq = g->nq;
    int nx2 = g->nx2;
    int ng11 = g->ng11;
    int d1 = g->d1;
    int d2 = g->d2;

    int i, j, q;
    int icopy = ng11;

    for(i=0; i<ng11; i++)
        for(j=0; j<nx2; j++)
            for(q=0; q<nq; q++)
                g->prim[d1*i+d2*j+q] = g->prim[d1*icopy+d2*j+q];
}

void bc_1R_outflow(struct grid *g, struct parList *par)
{
    int nq = g->nq;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int ng12 = g->ng12;
    int d1 = g->d1;
    int d2 = g->d2;

    int i, j, q;
    int icopy = nx1-ng12-1;

    for(i=nx1-ng12; i<nx1; i++)
        for(j=0; j<nx2; j++)
            for(q=0; q<nq; q++)
                g->prim[d1*i+d2*j+q] = g->prim[d1*icopy+d2*j+q];
}

void bc_2L_outflow(struct grid *g, struct parList *par)
{
    int nq = g->nq;
    int nx1 = g->nx1;
    int ng21 = g->ng21;
    int d1 = g->d1;
    int d2 = g->d2;

    int i, j, q;
    int jcopy = ng21;

    for(i=0; i<nx1; i++)
        for(j=0; j<ng21; j++)
            for(q=0; q<nq; q++)
                g->prim[d1*i+d2*j+q] = g->prim[d1*i+d2*jcopy+q];
}

void bc_2R_outflow(struct grid *g, struct parList *par)
{
    int nq = g->nq;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int ng22 = g->ng22;
    int d1 = g->d1;
    int d2 = g->d2;

    int i, j, q;
    int jcopy = nx2-ng22-1;

    for(i=0; i<nx1; i++)
        for(j=nx2-ng22; j<nx2; j++)
            for(q=0; q<nq; q++)
                g->prim[d1*i+d2*j+q] = g->prim[d1*i+d2*jcopy+q];
}
