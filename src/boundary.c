#include <stdio.h>
#include "boundary.h"
#include "par.h"

int set_boundary(struct parList *pars)
{
    int err = 0;
    int choice = pars->bc1Inner;

    if(choice == 0)
        bc_1L = &bc_1L_none;
    else if(choice == 1)
        bc_1L = &bc_1L_fixed;
    else if(choice == 2)
        bc_1L = &bc_1L_outflow;
    else if(choice == 3)
        bc_1L = &bc_1L_periodic;
    else
    {
        err++;
        printf("ERROR - Invalid choice for BC 1L: %d\n", choice);
    }

    choice = pars->bc1Outer;

    if(choice == 0)
        bc_1R = &bc_1R_none;
    else if(choice == 1)
        bc_1R = &bc_1R_fixed;
    else if(choice == 2)
        bc_1R = &bc_1R_outflow;
    else if(choice == 3)
        bc_1R = &bc_1R_periodic;
    else
    {
        err++;
        printf("ERROR - Invalid choice for BC 1R: %d\n", choice);
    }
    
    choice = pars->bc2Inner;

    if(choice == 0)
        bc_2L = &bc_2L_none;
    else if(choice == 1)
        bc_2L = &bc_2L_fixed;
    else if(choice == 2)
        bc_2L = &bc_2L_outflow;
    else if(choice == 3)
        bc_2L = &bc_2L_periodic;
    else
    {
        err++;
        printf("ERROR - Invalid choice for BC 2L: %d\n", choice);
    }
    
    choice = pars->bc2Outer;

    if(choice == 0)
        bc_2R = &bc_2R_none;
    else if(choice == 1)
        bc_2R = &bc_2R_fixed;
    else if(choice == 2)
        bc_2R = &bc_2R_outflow;
    else if(choice == 3)
        bc_2R = &bc_2R_periodic;
    else
    {
        err++;
        printf("ERROR - Invalid choice for BC 2R: %d\n", choice);
    }

    return err;
}
