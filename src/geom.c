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
        geom_CM2 = &geom_CM2_cart_xyz;
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
        geom_CM2 = &geom_CM2_cyl_rpz;
        geom_dA = &geom_dA_cyl_rpz;
        geom_dV = &geom_dV_cyl_rpz;
        geom_J = &geom_J_cyl_rpz;
        geom_J2 = &geom_J2_cyl_rpz;
        geom_gam = &geom_gam_cyl_rpz;
        geom_igam = &geom_igam_cyl_rpz;
        geom_dgam = &geom_dgam_cyl_rpz;
    }
    else if(choice == 2)
    {
        geom_CM = &geom_CM_sph_tpr;
        geom_CM2 = &geom_CM2_sph_tpr;
        geom_dA = &geom_dA_sph_tpr;
        geom_dV = &geom_dV_sph_tpr;
        geom_J = &geom_J_sph_tpr;
        geom_J2 = &geom_J2_sph_tpr;
        geom_gam = &geom_gam_sph_tpr;
        geom_igam = &geom_igam_sph_tpr;
        geom_dgam = &geom_dgam_sph_tpr;
    }
    else if(choice == 3)
    {
        geom_CM = &geom_CM_cyl_rpz_B;
        geom_CM2 = &geom_CM2_cyl_rpz_B;
        geom_dA = &geom_dA_cyl_rpz_B;
        geom_dV = &geom_dV_cyl_rpz_B;
        geom_J = &geom_J_cyl_rpz_B;
        geom_J2 = &geom_J2_cyl_rpz_B;
        geom_gam = &geom_gam_cyl_rpz_B;
        geom_igam = &geom_igam_cyl_rpz_B;
        geom_dgam = &geom_dgam_cyl_rpz_B;
    }
    else if(choice == 4)
    {
        geom_CM = &geom_CM_sph_rtp;
        geom_CM2 = &geom_CM2_sph_rtp;
        geom_dA = &geom_dA_sph_rtp;
        geom_dV = &geom_dV_sph_rtp;
        geom_J = &geom_J_sph_rtp;
        geom_J2 = &geom_J2_sph_rtp;
        geom_gam = &geom_gam_sph_rtp;
        geom_igam = &geom_igam_sph_rtp;
        geom_dgam = &geom_dgam_sph_rtp;
    }
    else if(choice == 5)
    {
        geom_CM = &geom_CM_sph_rtp_B;
        geom_CM2 = &geom_CM2_sph_rtp_B;
        geom_dA = &geom_dA_sph_rtp_B;
        geom_dV = &geom_dV_sph_rtp_B;
        geom_J = &geom_J_sph_rtp_B;
        geom_J2 = &geom_J2_sph_rtp_B;
        geom_gam = &geom_gam_sph_rtp_B;
        geom_igam = &geom_igam_sph_rtp_B;
        geom_dgam = &geom_dgam_sph_rtp_B;
    }
    else
    {
        err = 1;
        printf("ERROR - Invalid geometry choice: %d\n", choice);
    }

    return err;
}

