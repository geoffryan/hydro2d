#include <stdio.h>
#include "geom.h"
#include "par.h"

int set_geometry(struct parList *par)
{
    int err = 0;
    int choice = par->geom;

    if(choice == 0)
    {
        geom_CM = &geom_CM_flat_xyz;
        geom_dA = &geom_dA_flat_xyz;
        geom_dV = &geom_dV_flat_xyz;
        geom_gam = &geom_gam_flat_xyz;
        geom_igam = &geom_igam_flat_xyz;
        geom_dgam = &geom_dgam_flat_xyz;
    }
    else
    {
        err = 1;
        printf("ERROR - Invalid geometry choice: %d\n", choice);
    }

    return err;
}

