#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
#include "geom.h"
#include "hydro.h"
#include "par.h"
#include "timestep.h"

int set_timestep(struct parList *pars)
{
    int err = 0;
    int choice = pars->step;

    if(choice == 0)
        timestep = &step_fe;
    else if(choice == 1)
        timestep = &step_rk2_mp;
    else if(choice == 2)
        timestep = &step_rk2_tvd;
    else if(choice == 3)
        timestep = &step_rk3_tvd;
    else
    {
        printf("ERROR - Invalid Timestep choice: %d\n", choice);
        err = 1;
    }

    return err;
}

double get_dt(struct grid *g, struct parList *pars)
{
    int i,j;
    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int nq = g->nq;

    double dtmin = 1.0e100;

    for(i=0; i<nx1; i++)
        for(j=0; j<nx2; j++)
        {
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i+1], g->x2[j+1]};
            double x[2];
            geom_CM(xm,xp,x);
           
            double dx[2] = {xp[0]-xm[0], xp[1]-xm[1]};

            double dt = mindt(&(g->prim[nq*(nx1*i+j)]), x, dx, pars);

            dtmin = dt < dtmin ? dt : dtmin;
        }

    return pars->cfl * dtmin;
}

void step_fe(struct grid *g, double dt, struct parList *pars)
{
    substep(g, 1.0, 0.0, dt, pars);
    g->t += dt;
}

void step_rk2_mp(struct grid *g, double dt, struct parList *pars)
{
    // The Midpoint RK2 Method.

    copy_to_rk(g);
    substep(g, 0.5, 0.5, 0.5*dt, pars);
    g->t += 0.5*dt;

    update_cons(g, 0.0, 1.0);
    substep(g, 0.0, 1.0, dt, pars);
    g->t += 0.5*dt;
}

void step_rk2_tvd(struct grid *g, double dt, struct parList *pars)
{
    // The TVD RK2 Method by Gottlied & Shu

    double RK;

    RK = 1.0;
    copy_to_rk(g);
    substep(g, RK, 1.0-RK, RK*dt, pars);

    RK = 0.5;
    update_cons(g, RK, 1.0-RK);
    substep(g, RK, 1.0-RK, RK*dt, pars);
    g->t += dt;
}

void step_rk3_tvd(struct grid *g, double dt, struct parList *pars)
{
    // The TVD RK3 Method by Shu & Osher

    double RK;

    RK = 1.0;
    copy_to_rk(g);
    substep(g, RK, 1.0-RK, RK*dt, pars);

    RK = 0.25;
    update_cons(g, RK, 1.0-RK);
    substep(g, RK, 1.0-RK, RK*dt, pars);

    RK = 2.0/3.0;
    update_cons(g, RK, 1.0-RK);
    substep(g, RK, 1.0-RK, RK*dt, pars);

    g->t += dt;
}

void calc_cons(struct grid *g, struct parList *pars)
{
    int i, j;
    int nq = g->nq;
    int nx1 = g->nx1;
    int nx2 = g->nx2;

    for(i=0; i<nx1; i++)
        for(j=0; j<nx2; j++)
        {
            double xm[2] = {g->x1[i], g->x2[j]};
            double xp[2] = {g->x1[i+1], g->x2[j+1]};
            double x[2];
            geom_CM(xm,xp,x);
            double dV = geom_dV(xm,xp);

            prim2cons(&(g->prim[nq*(nx2*i+j)]), &(g->cons[nq*(nx2*i+j)]), 
                        x, dV, pars);
        }
}

void calc_prim(struct grid *g, struct parList *pars)
{
    int i, j;
    int nq = g->nq;
    int nx1 = g->nx1;
    int nx2 = g->nx2;

    double xm[2] = {0,0};
    double xp[2] = {0,0};
    double x[2] = {0,0};

    for(i=0; i<nx1; i++)
        for(j=0; j<nx2; j++)
        {
            xm[0] = g->x1[i];
            xm[1] = g->x2[j];
            xp[0] = g->x1[i+1];
            xp[1] = g->x2[j+1];
            geom_CM(xm,xp,x);
            double dV = geom_dV(xm,xp);

            cons2prim(&(g->prim[nq*(nx2*i+j)]), &(g->cons[nq*(nx2*i+j)]), 
                        x, dV, pars);
        }
}

