#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "geom.h"
#include "grid.h"
#include "hydro.h"
#include "par.h"

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
    g->ng11 = pars->bc1Inner == 0 ? 0 : pars->nghost;
    g->ng12 = pars->bc1Outer == 0 ? 0 : pars->nghost;
    g->ng21 = pars->bc2Inner == 0 ? 0 : pars->nghost;
    g->ng22 = pars->bc2Outer == 0 ? 0 : pars->nghost;
    g->nx1 = pars->nx1 + g->ng11 + g->ng12;
    g->nx2 = pars->nx2 + g->ng21 + g->ng22;
    g->x1min = pars->x1min;
    g->x1max = pars->x1max;
    g->x2min = pars->x2min;
    g->x2max = pars->x2max;
    g->t = pars->tmin;
    g->nc = pars->nc;
    g->nq = pars->nc + pars->np;
    g->PLM = pars->plm;

    g->d1 = g->nx2 * g->nq;
    g->d2 = g->nq;

    g->x1 = (double *) malloc((g->nx1+1) * sizeof(double));
    g->x2 = (double *) malloc((g->nx2+1) * sizeof(double));
    g->prim = (double *) malloc(g->nx1 * g->nx2 * g->nq * sizeof(double));
    g->cons = (double *) malloc(g->nx1 * g->nx2 * g->nq * sizeof(double));
    g->cons_rk = (double *) malloc(g->nx1 * g->nx2 * g->nq * sizeof(double));

    int i;
    double dx1 = (g->x1max - g->x1min) / g->nx1_int;
    double dx2 = (g->x2max - g->x2min) / g->nx2_int;
    for(i=0; i<g->nx1+1; i++)
        g->x1[i] = g->x1min + (i - g->ng11) * dx1;
    for(i=0; i<g->nx2+1; i++)
        g->x2[i] = g->x2min + (i - g->ng21) * dx2;
}

void free_grid(struct grid *g)
{
    free(g->x1);
    free(g->x2);
    free(g->prim);
    free(g->cons);
    free(g->cons_rk);
}

void interpolate_constant(struct grid *g, int i, int j, int dir,
                        double primL[], double primR[], struct parList *par)
{
    int q;
    int nq = g->nq;
    int d1 = g->d1;
    int d2 = g->d2;
    int iL, jL, iR, jR;

    if(dir == 0)
    {
        iL = i-1;
        iR = i;
        jL = j;
        jR = j;
    }
    else if(dir == 1)
    {
        iL = i;
        iR = i;
        jL = j-1;
        jR = j;
    }
    else
    {
        printf("ERROR - interpolate constant has bad dir=%d\n", dir);
        return;
    }

    for(q=0; q<nq; q++)
    {
        primL[q] = g->prim[d1*iL + d2*jL + q];
        primR[q] = g->prim[d1*iR + d2*jR + q];
    }
}

void interpolate_plm(struct grid *g, int i, int j, int dir,
                        double primL[], double primR[], struct parList *par)
{
    int q;
    int nq = g->nq;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int d1 = g->d1;
    int d2 = g->d2;
    int iLL, jLL, iL, jL, iR, jR, iRR, jRR, iRRR, jRRR;
    double plm = par->plm;

    if(dir == 0)
    {
        iLL = i-2;
        iL = i-1;
        iR = i;
        iRR = i+1;
        iRRR = i+2;
        jLL = j;
        jL = j;
        jR = j;
        jRR = j;
        jRRR = j;

        iLL = iLL >= 0 ? iLL : 0;
        iRRR = iRRR <= nx1 ? iRRR : nx1;
    }
    else if(dir == 1)
    {
        iLL = i;
        iL = i;
        iR = i;
        iRR = i;
        iRRR = i;
        jLL = j-2;
        jL = j-1;
        jR = j;
        jRR = j+1;
        jRRR = j+2;

        jLL = jLL >= 0 ? jLL : 0;
        jRRR = jRRR <= nx2 ? jRRR : nx2;
    }
    else
    {
        printf("ERROR - interpolate plm has bad dir=%d\n", dir);
        return;
    }

    double xfLL[2] = {g->x1[iLL], g->x2[jLL]};
    double xfL[2] = {g->x1[iL], g->x2[jL]};
    double xfC[2] = {g->x1[iR], g->x2[jR]};
    double xfR[2] = {g->x1[iRR], g->x2[jRR]};
    double xfRR[2] = {g->x1[iRRR], g->x2[jRRR]};
    double xLL[2], xL[2], xR[2], xRR[2];
    geom_CM(xfLL, xfL, xLL);
    geom_CM(xfL, xfC, xL);
    geom_CM(xfC, xfR, xR);
    geom_CM(xfR, xfRR, xRR);
    
    if((dir==0 && iLL==iL) || (dir==1 && jLL==jL))
        xLL[dir] = 2 * xfL[dir] - xL[dir];
    if((dir==0 && iRRR==iRR) || (dir==1 && jRRR==jRR))
        xRR[dir] = 2 * xfR[dir] - xR[dir];

    double idxLL = 1.0/(xL[dir]-xLL[dir]);
    double idxLC = 1.0/(xR[dir]-xLL[dir]);
    double idxC = 1.0/(xR[dir]-xL[dir]);
    double idxRC = 1.0/(xRR[dir]-xL[dir]);
    double idxRR = 1.0/(xRR[dir]-xR[dir]);

    double sL, sC, sR, gradL[nq], gradR[nq];
    int refl[nq];
    reflectInds(dir, refl, nq);

    if((dir==0 && iLL==iL) || (dir==1 && jLL==jL))
        for(q=0; q<nq; q++)
        {
            if(refl[q])
            {
                sL = (g->prim[d1*iR+d2*jR+q] + g->prim[d1*iL+d2*jL+q])*idxLL;
                sR = (g->prim[d1*iR+d2*jR+q] - g->prim[d1*iL+d2*jL+q])*idxC;
                sC = (g->prim[d1*iR+d2*jR+q] + g->prim[d1*iL+d2*jL+q])*idxLC;
            }
            else
            {
                sL = 0.0;
                sR = (g->prim[d1*iR+d2*jR+q] - g->prim[d1*iL+d2*jL+q])*idxC;
                sC = (g->prim[d1*iR+d2*jR+q] - g->prim[d1*iL+d2*jL+q])*idxLC;
            }
            gradL[q] = minmod(plm*sL, sC, plm*sR);
        }
    else
        for(q=0; q<nq; q++)
        {
            sL = (g->prim[d1*iL+d2*jL+q] - g->prim[d1*iLL+d2*jLL+q])*idxLL;
            sC = (g->prim[d1*iR+d2*jR+q] - g->prim[d1*iLL+d2*jLL+q])*idxLC;
            sR = (g->prim[d1*iR+d2*jR+q] - g->prim[d1*iL+d2*jL+q])*idxC;
            gradL[q] = minmod(plm*sL, sC, plm*sR);
        }

    if((dir==0 && iRRR==iRR) || (dir==1 && jRRR==jRR))
        for(q=0; q<nq; q++)
        {
            if(refl[q])
            {
                sL = (g->prim[d1*iR+d2*jR+q] - g->prim[d1*iL+d2*jL+q])*idxC;
                sR = (-g->prim[d1*iR+d2*jR+q] - g->prim[d1*iL+d2*jL+q])*idxRR;
                sC = (-g->prim[d1*iR+d2*jR+q] - g->prim[d1*iL+d2*jL+q])*idxRC;
            }
            else
            {
                sL = (g->prim[d1*iR+d2*jR+q] - g->prim[d1*iL+d2*jL+q])*idxC;
                sR = 0.0;
                sC = (g->prim[d1*iR+d2*jR+q] - g->prim[d1*iL+d2*jL+q])*idxRC;
            }
            gradR[q] = minmod(plm*sL, sC, plm*sR);
        }
    else
        for(q=0; q<nq; q++)
        {
            sL = (g->prim[d1*iR+d2*jR+q] - g->prim[d1*iL+d2*jL+q])*idxC;
            sC = (g->prim[d1*iRR+d2*jRR+q] - g->prim[d1*iL+d2*jL+q])*idxRC;
            sR = (g->prim[d1*iRR+d2*jRR+q] - g->prim[d1*iR+d2*jR+q])*idxRR;
            gradR[q] = minmod(plm*sL, sC, plm*sR);
        }

    for(q=0; q<nq; q++)
    {
        primL[q] = g->prim[d1*iL+d2*jL+q] + gradL[q]*(xfC[dir]-xL[dir]);
        primR[q] = g->prim[d1*iR+d2*jR+q] + gradR[q]*(xfC[dir]-xR[dir]);
    }
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
