#include <stdio.h>
#include "boundary.h"
#include "geom.h"
#include "grid.h"
#include "hydro.h"
#include "initial.h"
#include "io.h"
#include "par.h"
#include "riemann.h"
#include "timestep.h"

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

    err += set_boundary(&pars);
    err += set_geometry(&pars);
    err += set_hydro(&pars);
    err += set_initial(&pars);
    err += set_io(&pars);
    err += set_reconstruction(&pars);
    err += set_riemann(&pars);
    err += set_timestep(&pars);

    if(err)
    {
        printf("Error during setup: %d\n", err);
        return 0;
    }

    make_grid(&grid, &pars);
    initialize_grid(&grid, &pars);
    io_init(&io, &pars);

    io_out(&io, &grid, &pars);

    int i = 1;
    while(grid.t < pars.tmax)
    {
        double dt = get_dt(&grid, &pars);
        dt = grid.t+dt > pars.tmax ? pars.tmax-grid.t : dt;

        printf("t: %.6e dt: %.6e\n", grid.t, dt);

        timestep(&grid, dt, &pars);
        io_out(&io, &grid, &pars);

        i++;
    }

    free_grid(&grid);
    printf("CALCULATIONS CORRECT\n");

    return 0;
}
