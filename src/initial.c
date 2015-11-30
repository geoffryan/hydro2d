#include <stdlib.h>
#include <stdio.h>
#include "boundary.h"
#include "geom.h"
#include "grid.h"
#include "initial.h"
#include "par.h"
#include "timestep.h"

int set_initial(struct parList *par)
{
    int err = 0;
    int choice = par->init;

    if(choice == 0)
        initial_value = &initial_uniform;
    else
    {
        err++;
        printf("ERROR - Invalid choice for Initial Condition: %d\n", choice);
    }

    return err;
}

void initialize_grid(struct grid *g, struct parList *par)
{
    int nq = g->nq;
    int nx2 = g->nx2;
    int ng11 = g->ng11;
    int ng21 = g->ng21;
    int nx1_int = g->nx1_int;
    int nx2_int = g->nx2_int;

    int i, j;

    for(i=ng11; i < ng11+nx1_int; i++)
        for(j=ng21; j < ng21+nx2_int; j++)
        {
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i+1], g->x2[j+1]};
            double x[2];
            geom_CM(xm, xp, x);
            
            initial_value(&(g->prim[nq*(nx2*i+j)]), x, nq, par);
        }

    bc_1L(g, par);
    bc_1R(g, par);
    bc_2L(g, par);
    bc_2R(g, par);

    calc_cons(g, par);
}

