#include <stdio.h>
#include "boundary.h"
#include "geom.h"
#include "grid.h"
#include "par.h"
#include "initial.h"
#include "io.h"

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        printf("usage: hydro2d <parameter file>\n");
        return 0;
    }
    printf("\nWelcome to Hydro2D.\n\n");
    printf("Git Hash: %s\n", VERSION);

    int err = 0;
    struct parList pars = PAR_DEFAULT;
    struct grid grid = GRID_DEFAULT;
    struct io io = IO_DEFAULT;

    read_pars(&pars, argv[1]);

    err += set_reconstruction(&pars);
    err += set_geometry(&pars);
    err += set_io(&pars);
    err += set_boundary(&pars);

    if(err)
    {
        printf("Error during setup: %d\n", err);
        return 0;
    }

    make_grid(&grid, &pars);
    initialize_grid(&grid, &pars);

    io_out(&io, &grid, &pars);
    //print_pars(&pars, NULL);

    free_grid(&grid);

    return 0;
}
