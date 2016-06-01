#include <math.h>
#include "../geom.h"

void geom_CM_sph_tpr(double xm[], double xp[], double xc[])
{
    xc[0] = 0.5*(xm[0]+xp[0]);
    xc[1] = 0.5*(xm[1]+xp[1]);
}

void geom_CM2_sph_tpr(double xm[], double xp[], double xc[], int dir)
{
    xc[0] = 0.5*(xm[0]+xp[0]);
    xc[1] = 0.5*(xm[1]+xp[1]);
}

double geom_dA_sph_tpr(double xm[], double xp[], int dir)
{
    double dA;

    if(dir == 0)
        dA = sin(xm[0])*(xp[1]-xm[1]);
    else if(dir == 1)
        dA = xp[0]-xm[0];
    else
        dA = 0.0;

    return dA;
}

double geom_dV_sph_tpr(double xm[], double xp[])
{
    return (cos(xm[0])-cos(xp[0]))*(xp[1]-xm[1]);
}

double geom_J_sph_tpr(double x[])
{
    return sin(x[0]);
}

double geom_J2_sph_tpr(double x[], int dir)
{
    if(dir == 0)
        return sin(x[0]);
    else if(dir == 1)
        return 1.0;
    return 0.0;
}

void geom_gam_sph_tpr(double x[], double gam[3][3])
{
    double st = sin(x[0]);
    gam[0][0] = 1.0;
    gam[0][1] = 0.0;
    gam[0][2] = 0.0;
    gam[1][0] = 0.0;
    gam[1][1] = st*st;
    gam[1][2] = 0.0;
    gam[2][0] = 0.0;
    gam[2][1] = 0.0;
    gam[2][2] = 1.0;
}

void geom_igam_sph_tpr(double x[], double igam[3][3])
{
    double st = sin(x[0]);
    igam[0][0] = 1.0;
    igam[0][1] = 0.0;
    igam[0][2] = 0.0;
    igam[1][0] = 0.0;
    igam[1][1] = 1.0/(st*st);
    igam[1][2] = 0.0;
    igam[2][0] = 0.0;
    igam[2][1] = 0.0;
    igam[2][2] = 1.0;
}

void geom_dgam_sph_tpr(double x[], double dgam[2][3][3])
{
    int k,i,j;
    for(k=0; k<2; k++)
        for(i=0; i<3; i++)
            for(j=0; j<3; j++)
                dgam[k][i][j] = 0.0;
    dgam[0][1][1] = 2*sin(x[0])*cos(x[0]);
}
