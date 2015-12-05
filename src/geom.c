#include <stdio.h>
#include "geom.h"
#include "par.h"

int set_geometry(struct parList *par)
{
    int err = 0;
    int choice = par->geom;

    if(choice == 0)
    {
        geom_CM = &geom_CM_cart_xyz;
        geom_dA = &geom_dA_cart_xyz;
        geom_dV = &geom_dV_cart_xyz;
        geom_J = &geom_J_cart_xyz;
        geom_J2 = &geom_J2_cart_xyz;
        geom_gam = &geom_gam_cart_xyz;
        geom_igam = &geom_igam_cart_xyz;
        geom_dgam = &geom_dgam_cart_xyz;
    }
    else if(choice == 1)
    {
        geom_CM = &geom_CM_cyl_rpz;
        geom_dA = &geom_dA_cyl_rpz;
        geom_dV = &geom_dV_cyl_rpz;
        geom_J = &geom_J_cyl_rpz;
        geom_J2 = &geom_J2_cyl_rpz;
        geom_gam = &geom_gam_cyl_rpz;
        geom_igam = &geom_igam_cyl_rpz;
        geom_dgam = &geom_dgam_cyl_rpz;
    }
    else
    {
        err = 1;
        printf("ERROR - Invalid geometry choice: %d\n", choice);
    }

    return err;
}

