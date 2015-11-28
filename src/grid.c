#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "par.h"
#include "grid.h"

//Local Functions
double minmod(double a, double b, double c);

int set_reconstruction(struct parList *pars)
{
    int err = 0;
    int choice = pars->recon;

    if(choice == 0)
        reconstruction = &interpolate_constant;
    else if(choice == 1)
        reconstruction = &interpolate_plm;
    else
    {
        err = 1;
        printf("ERROR - Invalid Reconstruction choice: %d\n", choice);
    }

    return err;
}

void make_grid(struct grid *g, struct parList *pars)
{
    //Set grid variables from parameter list and allocate
    //space for data arrays.
    g->nx1_int = pars->nx1;
    g->nx2_int = pars->nx2;
    g->nx1 = pars->nx1 + 2*pars->nghost;
    g->nx2 = pars->nx2 + 2*pars->nghost;
    g->ng = pars->nghost;
    g->x1min = pars->x1min;
    g->x1max = pars->x1max;
    g->x2min = pars->x2min;
    g->x2max = pars->x2max;
    g->t = pars->tmin;
    g->nc = pars->nc;
    g->nq = pars->nc + pars->np;
    g->PLM = pars->plm;

    g->x1 = (double *) malloc((g->nx1+1) * sizeof(double));
    g->x2 = (double *) malloc((g->nx2+1) * sizeof(double));
    g->prim = (double *) malloc(g->nx1 * g->nx2 * g->nq * sizeof(double));
    g->cons = (double *) malloc(g->nx1 * g->nx2 * g->nq * sizeof(double));
    g->cons_rk = (double *) malloc(g->nx1 * g->nx2 * g->nq * sizeof(double));
    g->grad = (double *) malloc(g->nx1 * g->nx2 * g->nq * 2 * sizeof(double));

    int i;
    double dx1 = (g->x1max - g->x1min) / g->nx1_int;
    double dx2 = (g->x2max - g->x2min) / g->nx2_int;
    for(i=0; i<g->nx1+1; i++)
        g->x1[i] = g->x1min + (i - g->ng) * dx1;
    for(i=0; i<g->nx2+1; i++)
        g->x2[i] = g->x2min + (i - g->ng) * dx2;
}

void free_grid(struct grid *g)
{
    free(g->x1);
    free(g->x2);
    free(g->prim);
    free(g->cons);
    free(g->cons_rk);
    free(g->grad);
}

void interpolate_constant(struct grid *g)
{
    int i,q;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int nq = g->nq;

    for(i=0; i<nx1*nx2; i++)
        for(q=0; q<nq; q++)
        {
            g->grad[2*nq*i+q+0] = 0.0;
            g->grad[2*nq*i+q+1] = 0.0;
        }
}

void interpolate_plm(struct grid *g)
{
    //TODO: Write this.
    /*
    int i,j;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int nq = g->nq;
    double PLM = g->PLM;

    for(j=0; j<nq; j++)
    {
        double dR = 2*(g->prim[1*nq+j] - g->prim[0*nq+j]) / (g->x[2]-g->x[0]);
        g->grad[j] = dR;
    }
    for(i=1; i<nx-1; i++)
    {
        double idxL = 2.0 / (g->x[i+1] - g->x[i-1]);
        double idxC = 2.0 / (g->x[i+2] + g->x[i+1] - g->x[i] - g->x[i-1]);
        double idxR = 2.0 / (g->x[i+2] - g->x[i]);
        int id = nq*i;
        
        for(j=0; j<nq; j++)
        {
            double SL = (g->prim[id   +j] - g->prim[id-nq+j]) * idxL;
            double SC = (g->prim[id+nq+j] - g->prim[id-nq+j]) * idxC;
            double SR = (g->prim[id+nq+j] - g->prim[id   +j]) * idxR;
            
            g->grad[id+j] = minmod(PLM*SL, SC, PLM*SR);
        }
    }
    for(j=0; j<nq; j++)
    {
        double dL = 2*(g->prim[(nx-1)*nq+j] - g->prim[(nx-2)*nq+j])
                        / (g->x[nx]-g->x[nx-2]);
        g->grad[j] = dL;
    }
    */
}

void copy_to_rk(struct grid *g)
{
    // Copy cons into cons_rk.

    int i;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int nq = g->nq;

    for(i=0; i<nx1*nx2*nq; i++)
        g->cons_rk[i] = g->cons[i];
}

void update_cons(struct grid *g, double fac1, double fac2)
{
    // Update cons with fac1*cons & fac2*cons_rk
    
    int i;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int nq = g->nq;

    for(i=0; i<nx1*nx2*nq; i++)
        g->cons[i] = fac1 * g->cons[i] + fac2 * g->cons_rk[i];
}

void update_cons_rk(struct grid *g, double fac1, double fac2)
{
    // Update cons_rk with fac1*cons & fac2*cons_rk

    int i;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int nq = g->nq;

    for(i=0; i<nx1*nx2*nq; i++)
        g->cons_rk[i] = fac1 * g->cons[i] + fac2 * g->cons_rk[i];
}

//Local definitions

double minmod(double a, double b, double c)
{
    double m;

    if(a*b < 0)
        m = 0;
    else if(fabs(a) < fabs(b))
        m = a;
    else
        m = b;

    if(m*c < 0)
        m = 0;
    else if(fabs(c) < fabs(m))
        m = c;

    return m;
}
