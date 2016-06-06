#include <math.h>
#include "../geom.h"

void geom_CM_sph_rtp(double xm[], double xp[], double xc[])
{
    xc[0] = 0.5*(xm[0]+xp[0]);
    xc[1] = 0.5*(xm[1]+xp[1]);
}

void geom_CM2_sph_rtp(double xm[], double xp[], double xc[], int dir)
{
    xc[0] = 0.5*(xm[0]+xp[0]);
    xc[1] = 0.5*(xm[1]+xp[1]);
}

double geom_dA_sph_rtp(double xm[], double xp[], int dir)
{
    double dA;

    if(dir == 0)
        dA = xm[0]*(cos(xm[1])-cos(xp[1]));
    else if(dir == 1)
        dA = 0.5*(xp[0]+xm[0])*(xp[0]-xm[0])*(cos(xm[1])-cos(xp[1]));
    else
        dA = 0.0;

    return dA;
}

double geom_dV_sph_rtp(double xm[], double xp[])
{
    return (xp[0]+xm[0])*(xp[0]-xm[0])*(cos(xm[1])-cos(xp[1]))/3.0;
}

double geom_J_sph_rtp(double x[])
{
    return x[0]*x[0]*sin(x[1]);
}

double geom_J2_sph_rtp(double x[], int dir)
{
    if(dir == 0)
        return x[0]*x[0]*sin(x[1]);
    else if(dir == 1)
        return x[0]*sin(x[1]);
    return 0.0;
}

void geom_gam_sph_rtp(double x[], double gam[3][3])
{
    double st = sin(x[0]);
    gam[0][0] = 1.0;
    gam[0][1] = 0.0;
    gam[0][2] = 0.0;
    gam[1][0] = 0.0;
    gam[1][1] = x[0]*x[0];
    gam[1][2] = 0.0;
    gam[2][0] = 0.0;
    gam[2][1] = 0.0;
    gam[2][2] = x[0]*x[0]*st*st;
}

void geom_igam_sph_rtp(double x[], double igam[3][3])
{
    double st = sin(x[0]);
    igam[0][0] = 1.0;
    igam[0][1] = 0.0;
    igam[0][2] = 0.0;
    igam[1][0] = 0.0;
    igam[1][1] = 1.0/(x[0]*x[0]);
    igam[1][2] = 0.0;
    igam[2][0] = 0.0;
    igam[2][1] = 0.0;
    igam[2][2] = 1.0/(x[0]*x[0]*st*st);
}

void geom_dgam_sph_rtp(double x[], double dgam[2][3][3])
{
    double st = sin(x[0]);
    double ct = cos(x[0]);
    int k,i,j;
    for(k=0; k<2; k++)
        for(i=0; i<3; i++)
            for(j=0; j<3; j++)
                dgam[k][i][j] = 0.0;

    dgam[0][1][1] = 2*x[0];
    dgam[0][2][2] = 2*x[0]*st*st;
    dgam[1][2][2] = 2*x[0]*x[0]*st*ct;
}
