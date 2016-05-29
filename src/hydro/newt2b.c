#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../par.h"
#include "../geom.h"
#include "../hydro.h"
#include "newt2b.h"

void prim2cons_newt2b(double *prim, double *cons, double x[2], double dV,
                            struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double lx1 = prim[LX1];
    double lx2 = prim[LX2];
    double adind = pars->gammalaw;

    double igam[3][3];
    geom_igam(x, igam);

    double v2 = igam[0][0]*lx1*lx1 + igam[1][1]*lx2*lx2 + 2*igam[0][1]*lx1*lx2;

    cons[DDD] = rho * dV;
    cons[TAU] = (0.5*rho*v2 + P/(adind-1)) * dV;
    cons[SX1] = rho * lx1 * dV;
    cons[SX2] = rho * lx2 * dV;
}

void cons2prim_newt2b(double *cons, double *prim, double x[2], double dV,
                            struct parList *pars)
{
    double rho = cons[DDD] / dV;
    double en = cons[TAU] / dV;
    double mx1 = cons[SX1] / dV;
    double mx2 = cons[SX2] / dV;
    double adind = pars->gammalaw;
    
    double igam[3][3];
    geom_igam(x, igam);

    double m2 = igam[0][0]*mx1*mx1 + igam[1][1]*mx2*mx2 + 2*igam[0][1]*mx1*mx2;

    prim[RHO] = rho;
    prim[PPP] = (adind-1.0)*(en - 0.5*m2/rho);
    prim[LX1] = mx1 / rho;
    prim[LX2] = mx2 / rho;
}

void flux_newt2b(double *prim, double *F, double x[2], int dir,
                struct parList *pars)
{
    double rho = prim[RHO];
    double P = prim[PPP];
    double lx1 = prim[LX1];
    double lx2 = prim[LX2];
    double adind = pars->gammalaw;

    double igam[3][3];
    geom_igam(x, igam);

    double v2 = igam[0][0]*lx1*lx1 + igam[1][1]*lx2*lx2 + 2*igam[0][1]*lx1*lx2;
    double vx1 = igam[0][0]*lx1 + igam[0][1]*lx2;
    double vx2 = igam[1][0]*lx1 + igam[1][1]*lx2;

    if(dir == 0)
    {
        F[DDD] = rho * vx1;
        F[TAU] = (0.5*rho*v2 + adind/(adind-1.0)*P) * vx1;
        F[SX1] = rho*lx1*vx1 + P;
        F[SX2] = rho*lx2*vx1;
    }
    else if(dir == 1)
    {
        F[DDD] = rho * vx2;
        F[TAU] = (0.5*rho*v2 + adind/(adind-1.0)*P) * vx2;
        F[SX1] = rho*lx1*vx2;
        F[SX2] = rho*lx2*vx2 + P;
    }
    else
    {
        printf("ERROR in newt2 Flux. Bad dir=%d\n", dir);
    }
}

void add_source_newt2b(double *prim, double *cons, double xm[2], double xp[2],
                        double dt, struct parList *pars)
{
    double x[2];
    geom_CM(xm, xp, x);

    double igam[3][3];
    double dgam[2][3][3];
    geom_igam(x, igam);
    geom_dgam(x, dgam);

    double rho = prim[RHO];
    double lx[2] = {prim[LX1], prim[LX2]};
    double vx[2] = {0.0, 0.0};
    double P = prim[PPP];

    int i, j;
    double sx1 = 0;
    double sx2 = 0;

    for(i=0; i<2; i++)
        for(j=0; j<2; j++)
            vx[i] += igam[i][j]*lx[j];

    for(i=0; i<2; i++)
        for(j=0; j<2; j++)
        {
            sx1 += 0.5 * (rho*vx[i]*vx[j] + igam[i][j]*P) * dgam[0][i][j];
            sx2 += 0.5 * (rho*vx[i]*vx[j] + igam[i][j]*P) * dgam[1][i][j];
        }

    cons[SX1] += sx1 * dVdt;
    cons[SX2] += sx2 * dVdt;
}

void wave_speeds_newt2b(double *primA, double *primB, double *sL, 
                        double *sR, double *sC, double x[2], int dir,
                        struct parList *pars)
{
    double hh, igam[3][3];
    geom_igam(x, igam);
    hh = sqrt(igam[dir][dir]);

    double rhoA = primA[RHO];
    double PA = primA[PPP];
    double lx1A = primA[LX1];
    double lx2A = primA[LX2];
    double rhoB = primB[RHO];
    double PB = primB[PPP];
    double lx1B = primB[LX1];
    double lx2B = primB[LX2];
    double adind = pars->gammalaw;
    
    double csA = sqrt(adind*PA/rhoA) * hh;
    double csB = sqrt(adind*PB/rhoB) * hh;
    
    double vx1A = igam[0][0]*lx1A + igam[0][1]*lx2A;
    double vx2A = igam[1][0]*lx1A + igam[1][1]*lx2A;
    double vx1B = igam[0][0]*lx1B + igam[0][1]*lx2B;
    double vx2B = igam[1][0]*lx1B + igam[1][1]*lx2B;

    if(dir == 0)
    {
        *sL = vx1A - csA;
        if(vx1B - csB < *sL)
            *sL = vx1B - csB;

        *sR = vx1A + csA;
        if(vx1B + csB > *sR)
            *sR = vx1B + csB;

        //*sC = (P2-P1 - rho1*v11*cs1 - rho2*v21*cs2) / (-rho1*cs1-rho2*cs2);
    }
    else if(dir == 1)
    {
        *sL = vx2A - csA;
        if(vx2B - csB < *sL)
            *sL = vx2B - csB;

        *sR = vx2A + csA;
        if(vx2B + csB > *sR)
            *sR = vx2B + csB;

        //*sC = (P2-P1 - rho1*v12*cs1 - rho2*v22*cs2) / (-rho1*cs1-rho2*cs2);
    }
    else
    {
        printf("ERROR in wavespeeds newt2: bad dir=%d\n", dir);
    }
}

double mindt_newt2b(double *prim, double x[2], double dx[2], 
                    struct parList *pars)
{
    double hh1, hh2, igam[3][3];
    geom_igam(x, igam);
    hh1 = sqrt(igam[0][0]);
    hh2 = sqrt(igam[1][1]);

    double rho = prim[RHO];
    double P = prim[PPP];
    double lx1 = prim[LX1];
    double lx2 = prim[LX2];
    double adind = pars->gammalaw;

    double vx1 = igam[0][0]*lx1 + igam[0][1]*lx2;
    double vx2 = igam[1][0]*lx1 + igam[1][1]*lx2;
    double cs = sqrt(adind*P/rho);
    
    double s1l = fabs(vx1 - cs*hh1);
    double s1r = fabs(vx1 + cs*hh1);
    double s2l = fabs(vx2 - cs*hh2);
    double s2r = fabs(vx2 + cs*hh2);

    double s1 = s1l > s1r ? s1l : s1r;
    double s2 = s2l > s2r ? s2l : s2r;

    double dt1 = dx[0] / s1;
    double dt2 = dx[1] / s2;

    double dt = dt1>dt2 ? dt2 : dt1;

    return dt;
}

void reflectInds_newt2b(int dir, int *inds, int nq)
{
    int i;
    for(i=0; i<nq; i++)
        inds[i] = 0;

    if(dir==0)
        inds[LX1] = 1;
    else if(dir==1)
        inds[LX2] = 1;
}

void Ustar_newt2b(double *prim, double *Us, double sK, double sC, double x[2],
                        struct parList *pars)
{
    /*
    double rho = prim[RHO];
    double P = prim[PPP];
    double vx1 = prim[VX1];
    double vx2 = prim[VX2];
    double adind = pars->gammalaw;

    double gam[3][3];
    geom_gam(x, gam);

    double v2 = gam[0][0]*vx1*vx1 + gam[1][1]*vx2*vx2 + 2*gam[0][1]*vx1*vx2;
    //double l1 = gam[0][0]*vx1 + gam[0][1]*vx2;
    //double l2 = gam[1][0]*vx1 + gam[1][1]*vx2;
    
    double rhostar = rho * (sK-vx1) / (sK - sC);
    double rhoestar = P/(adind-1.0) * (sK-vx1) / (sK - sC);
    double Pstar = P * (sC-vx1) / (sK - sC);

    Us[RHO] = rhostar;
    Us[TAU] = 0.5*rhostar*v2 + rhoestar + rhostar*sC*(sC-vx1) + Pstar;
    Us[SX1] = rhostar * sC;
    Us[SX2] = rhostar * vx2;
    */
}
