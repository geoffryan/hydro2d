#include "../geom.h"

void geom_CM_cart_xyz(double xm[], double xp[], double xc[])
{
    xc[0] = 0.5*(xm[0]+xp[0]);
    xc[1] = 0.5*(xm[1]+xp[1]);
}

void geom_CM2_cart_xyz(double xm[], double xp[], double xc[], int dir)
{
    xc[0] = 0.5*(xm[0]+xp[0]);
    xc[1] = 0.5*(xm[1]+xp[1]);
}

double geom_dA_cart_xyz(double xm[], double xp[], int dir)
{
    double dA;

    if(dir == 0)
        dA = xp[1]-xm[1];
    else if(dir == 1)
        dA = xp[0]-xm[0];
    else
        dA = 0.0;

    return dA;
}

double geom_dV_cart_xyz(double xm[], double xp[])
{
    return (xp[0]-xm[0])*(xp[1]-xm[1]);
}

double geom_J_cart_xyz(double x[])
{
    return 1.0;
}

double geom_J2_cart_xyz(double x[], int dir)
{
    return 1.0;
}

void geom_gam_cart_xyz(double x[], double gam[3][3])
{
    gam[0][0] = 1.0;
    gam[0][1] = 0.0;
    gam[0][2] = 0.0;
    gam[1][0] = 0.0;
    gam[1][1] = 1.0;
    gam[1][2] = 0.0;
    gam[2][0] = 0.0;
    gam[2][1] = 0.0;
    gam[2][2] = 1.0;
}

void geom_igam_cart_xyz(double x[], double igam[3][3])
{
    igam[0][0] = 1.0;
    igam[0][1] = 0.0;
    igam[0][2] = 0.0;
    igam[1][0] = 0.0;
    igam[1][1] = 1.0;
    igam[1][2] = 0.0;
    igam[2][0] = 0.0;
    igam[2][1] = 0.0;
    igam[2][2] = 1.0;
}

void geom_dgam_cart_xyz(double x[], double dgam[2][3][3])
{
    int k,i,j;
    for(k=0; k<2; k++)
        for(i=0; i<3; i++)
            for(j=0; j<3; j++)
                dgam[k][i][j] = 0.0;
}
