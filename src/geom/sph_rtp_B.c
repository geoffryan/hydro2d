#include <math.h>
#include <stdio.h>
#include "../geom.h"

void geom_CM_sph_rtp_B(double xm[], double xp[], double xc[])
{
    double r1 = xm[0];
    double r2 = xp[0];
    double th1 = xm[1];
    double th2 = xp[1];
    double ct1 = cos(th1);
    double ct2 = cos(th2);
    
    //cos(th1) - cos(th2)
    double dcost = 2*sin(0.5*(th1+th2))*sin(0.5*(th2-th1));
    //sin(th2) - sin(th1)
    double dsint = 2*cos(0.5*(th1+th2))*sin(0.5*(th2-th1));

    //xc[0] = 0.5*(r1+r2);
    xc[0] = 3.0*(r1*r1*r1+r1*r1*r2+r1*r2*r2+r2*r2*r2)
                / (4.0*(r1*r1+r1*r2+r2*r2));
    xc[1] = (th1*ct1 - th2*ct2 + dsint) / dcost;
}

void geom_CM2_sph_rtp_B(double xm[], double xp[], double xc[], int dir)
{
    double r1 = xm[0];
    double r2 = xp[0];
    double th1 = xm[1];
    double th2 = xp[1];
    double ct1 = cos(th1);
    double ct2 = cos(th2);
    
    //cos(th1) - cos(th2)
    double dcost = 2*sin(0.5*(th1+th2))*sin(0.5*(th2-th1));
    //sin(th2) - sin(th1)
    double dsint = 2*cos(0.5*(th1+th2))*sin(0.5*(th2-th1));

    //xc[0] = 0.5*(r1+r2);
    xc[0] = 3.0*(r1*r1*r1+r1*r1*r2+r1*r2*r2+r2*r2*r2)
                / (4.0*(r1*r1+r1*r2+r2*r2));
    if(dir == 0)
        xc[1] = (th1*ct1 - th2*ct2 + dsint) / dcost;
    else if(dir == 1)
        xc[1] = th1;
}

double geom_dA_sph_rtp_B(double xm[], double xp[], int dir)
{
    double dA;
    double r1 = xm[0];
    double r2 = xp[0];
    double th1 = xm[1];
    double th2 = xp[1];
    double st1 = sin(th1);

    double dcost = 2*sin(0.5*(th1+th2))*sin(0.5*(th2-th1));

    if(dir == 0)
        dA = r1*r1*dcost;
    else if(dir == 1)
        dA = (r1*r1+r1*r2+r2*r2)*(r2-r1)*st1/3.0;
    else
        dA = 0.0;

    return dA;
}

double geom_dV_sph_rtp_B(double xm[], double xp[])
{
    double r1 = xm[0];
    double r2 = xp[0];
    double th1 = xm[1];
    double th2 = xp[1];

    double dcost = 2*sin(0.5*(th1+th2))*sin(0.5*(th2-th1));

    return (r1*r1+r1*r2+r2*r2)*(r2-r1) * dcost / 3.0;
}

double geom_J_sph_rtp_B(double x[])
{
    return 1.0;
}

double geom_J2_sph_rtp_B(double x[], int dir)
{
    return 1.0;
}

void geom_gam_sph_rtp_B(double x[], double gam[3][3])
{
    double r = x[0];
    double st = sin(x[1]);
    gam[0][0] = 1.0;
    gam[0][1] = 0.0;
    gam[0][2] = 0.0;
    gam[1][0] = 0.0;
    gam[1][1] = r*r;
    gam[1][2] = 0.0;
    gam[2][0] = 0.0;
    gam[2][1] = 0.0;
    gam[2][2] = r*r*st*st;
}

void geom_igam_sph_rtp_B(double x[], double igam[3][3])
{
    double r = x[0];
    double st = sin(x[1]);
    igam[0][0] = 1.0;
    igam[0][1] = 0.0;
    igam[0][2] = 0.0;
    igam[1][0] = 0.0;
    igam[1][1] = 1.0/(r*r);
    igam[1][2] = 0.0;
    igam[2][0] = 0.0;
    igam[2][1] = 0.0;
    igam[2][2] = 1.0/(r*r*st*st);
}

void geom_dgam_sph_rtp_B(double x[], double dgam[2][3][3])
{
    double r = x[0];
    double st = sin(x[1]);
    double ct = cos(x[1]);
    int k,i,j;
    for(k=0; k<2; k++)
        for(i=0; i<3; i++)
            for(j=0; j<3; j++)
                dgam[k][i][j] = 0.0;

    dgam[0][1][1] = 2*r;
    dgam[0][2][2] = 2*r*st*st;
    dgam[1][2][2] = 2*r*r*st*ct;
}
