#include <stdlib.h>
#include <stdio.h>
#include "par.h"
#include "hydro.h"

int set_hydro(struct parList *pars)
{
    int err = 0;
    int choice = pars->hydro;

    if(choice == 0)
    {
        prim2cons = &prim2cons_newt2;
        cons2prim = &cons2prim_newt2;
        flux = &flux_newt2;
        add_source = &add_source_newt2;
        wave_speeds = &wave_speeds_newt2;
        mindt = &mindt_newt2;
        Ustar = &Ustar_newt2;
    }
    else
    {
        err++;
        printf("ERROR - Invalid choice for hydro: %d\n", choice);
    }

    return err;
}
