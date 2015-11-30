#include <stdlib.h>
#include "../boundary.h"
#include "../grid.h"
#include "../par.h"

//Periodic boundaries

void bc_1L_periodic(struct grid *g, struct parList *par)
{
    int nq = g->nq;
    int nx2 = g->nx2;
    int ng11 = g->ng11;

    int i, j, q;
    int iwrap = g->nx1_int;

    for(i=0; i<ng11; i++)
        for(j=0; j<nx2; j++)
            for(q=0; q<nq; q++)
                g->prim[nq*(nx2*i+j)+q] = g->prim[nq*(nx2*(iwrap+i)+j)+q];
}

void bc_1R_periodic(struct grid *g, struct parList *par)
{
    int nq = g->nq;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int ng12 = g->ng12;

    int i, j, q;
    int iwrap = g->ng11 - (nx1-ng12);

    for(i=nx1-ng12; i<nx1; i++)
        for(j=0; j<nx2; j++)
            for(q=0; q<nq; q++)
                g->prim[nq*(nx2*i+j)+q] = g->prim[nq*(nx2*(iwrap+i)+j)+q];
}

void bc_2L_periodic(struct grid *g, struct parList *par)
{
    int nq = g->nq;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int ng21 = g->ng21;

    int i, j, q;
    int jwrap = g->nx2_int;

    for(i=0; i<nx1; i++)
        for(j=0; j<ng21; j++)
            for(q=0; q<nq; q++)
                g->prim[nq*(nx2*i+j)+q] = g->prim[nq*(nx2*i+jwrap+j)+q];
}

void bc_2R_periodic(struct grid *g, struct parList *par)
{
    int nq = g->nq;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int ng22 = g->ng22;

    int i, j, q;
    int jwrap = g->ng21 - (nx2-ng22);

    for(i=0; i<nx1; i++)
        for(j=nx2-ng22; j<nx2; j++)
            for(q=0; q<nq; q++)
                g->prim[nq*(nx2*i+j)+q] = g->prim[nq*(nx2*i+jwrap+j)+q];
}
