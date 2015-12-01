#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include "hdf5.h"
#include "grid.h"
#include "par.h"
#include "io.h"

//Local Functions
void io_print_header_ascii(char filename[]);
void io_print_grid_ascii(struct grid *g, char filename[]);
void io_print_header_hdf5(hid_t id);
void io_print_pars_hdf5(struct parList *par, hid_t id);
void io_print_grid_hdf5(struct grid *g, struct parList *par, hid_t id);

//Utility Functions
void io_write_attr_hdf5(hid_t id, char name[], int type, void *val, int len);
void io_write_dset_hdf5(hid_t id, char name[], int type, void *dat, int rank,
                        hsize_t dims[]);

int set_io(struct parList *par)
{
    int err = 0;
    int choice = par->io;

    if(choice == 0)
        io_out = &io_out_ascii;
    else if(choice == 1)
        io_out = &io_out_hdf5;
    else
    {
        err = 1;
        printf("ERROR - Invalid IO choice: %d\n", choice);
    }

    return err;
}

void io_init(struct io *iom, struct parList *par)
{
    iom->nTot = par->nChkpt;
    iom->current = 0;
    iom->tnext = par->tmin;
    iom->tmin = par->tmin;
    iom->tmax = par->tmax;
}

void io_out_ascii(struct io *iom, struct grid *g, struct parList *par)
{
    if(g->t >= iom->tnext && iom->current <= iom->nTot)
    {
        char filename[128];

        sprintf(filename, "checkpoint_%05d.txt", iom->current);

        io_print_header_ascii(filename);
        print_pars(par, filename);
        io_print_grid_ascii(g, filename);

        iom->current += 1;
        iom->tnext = iom->current * (iom->tmax-iom->tmin)/(iom->nTot)
                        + iom->tmin;
        if(iom->current == iom->nTot)
            iom->tnext = iom->tmax;
    }
}

void io_out_hdf5(struct io *iom, struct grid *g, struct parList *par)
{
    if(g->t >= iom->tnext && iom->current <= iom->nTot)
    {
        char filename[128];

        sprintf(filename, "checkpoint_%05d.h5", iom->current);

        hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, 
                                H5P_DEFAULT);
        hid_t group = H5Gcreate(file, "/data", H5P_DEFAULT, H5P_DEFAULT, 
                                H5P_DEFAULT);

        io_print_header_hdf5(group);
        io_print_pars_hdf5(par, group);
        io_print_grid_hdf5(g, par, group);

        H5Gclose(group);
        H5Fclose(file);

        iom->current += 1;
        iom->tnext = iom->current * (iom->tmax-iom->tmin)/(iom->nTot)
                        + iom->tmin;
        if(iom->current == iom->nTot)
            iom->tnext = iom->tmax;
    }

}

void io_print_header_ascii(char filename[])
{
    FILE *f;
    if(filename == NULL)
        f = stdout;
    else
        f = fopen(filename, "w");

    fprintf(f, "### Hydro2D Checkpoint ###\n");
    fprintf(f, "GitHash %s\n", VERSION);
    
    if(filename != NULL)
        fclose(f);
}

void io_print_grid_ascii(struct grid *g, char filename[])
{
    FILE *f;
    if(filename == NULL)
        f = stdout;
    else
        f = fopen(filename, "a");

    int i,j,q;

    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int nq = g->nq;

    fprintf(f, "### Data ###\n");

    fprintf(f, "t %.12g\n", g->t);    
    fprintf(f, "x1");
    for(i=0; i<nx1+1; i++)
        fprintf(f, " %.12g", g->x1[i]);
    fprintf(f, "\n");

    fprintf(f, "x2");
    for(i=0; i<nx2+1; i++)
        fprintf(f, " %.12g", g->x2[i]);
    fprintf(f, "\n");

    for(i=0; i<nx1; i++)
        for(j=0; j < nx2; j++)
        {
            fprintf(f, "%d %d", i, j);
            for(q=0; q<nq; q++)
                fprintf(f, " %.12g", g->prim[nq*(nx2*i+j)+q]);
            for(q=0; q<nq; q++)
                fprintf(f, " %.12g", g->cons[nq*(nx2*i+j)+q]);
            fprintf(f, "\n");
        }

    if(filename != NULL)
        fclose(f);
}

void io_print_header_hdf5(hid_t id)
{
    io_write_attr_hdf5(id, "GitHash", VAR_STR, VERSION, sizeof(VERSION));
}

void io_print_pars_hdf5(struct parList *par, hid_t id)
{
    io_write_attr_hdf5(id, "Hydro", VAR_INT, &(par->hydro), 0);
    io_write_attr_hdf5(id, "Geom", VAR_INT, &(par->geom), 0);
    io_write_attr_hdf5(id, "EOS",   VAR_INT, &(par->eos), 0);
    io_write_attr_hdf5(id, "Cool",  VAR_INT, &(par->cool), 0);
    io_write_attr_hdf5(id, "Recon",  VAR_INT, &(par->recon), 0);
    io_write_attr_hdf5(id, "Riemann",  VAR_INT, &(par->riemann), 0);
    io_write_attr_hdf5(id, "Timestep",  VAR_INT, &(par->step), 0);
    io_write_attr_hdf5(id, "BCX1Inner",  VAR_INT, &(par->bc1Inner), 0);
    io_write_attr_hdf5(id, "BCX1Outer",  VAR_INT, &(par->bc1Outer), 0);
    io_write_attr_hdf5(id, "BCX2Inner",  VAR_INT, &(par->bc2Inner), 0);
    io_write_attr_hdf5(id, "BCX2Outer",  VAR_INT, &(par->bc2Outer), 0);
    io_write_attr_hdf5(id, "Nx1",     VAR_INT, &(par->nx1), 0);
    io_write_attr_hdf5(id, "Nx2",     VAR_INT, &(par->nx2), 0);
    io_write_attr_hdf5(id, "Nghost",     VAR_INT, &(par->nghost), 0);
    io_write_attr_hdf5(id, "Ncons",     VAR_INT, &(par->nc), 0);
    io_write_attr_hdf5(id, "Npass",     VAR_INT, &(par->np), 0);
    io_write_attr_hdf5(id, "Tmin",     VAR_DBL, &(par->tmin), 0);
    io_write_attr_hdf5(id, "Tmax",     VAR_DBL, &(par->tmax), 0);
    io_write_attr_hdf5(id, "X1min",     VAR_DBL, &(par->x1min), 0);
    io_write_attr_hdf5(id, "X1max",     VAR_DBL, &(par->x1max), 0);
    io_write_attr_hdf5(id, "X2min",     VAR_DBL, &(par->x2min), 0);
    io_write_attr_hdf5(id, "X2max",     VAR_DBL, &(par->x2max), 0);
    io_write_attr_hdf5(id, "PLM",     VAR_DBL, &(par->plm), 0);
    io_write_attr_hdf5(id, "CFL",     VAR_DBL, &(par->cfl), 0);
    io_write_attr_hdf5(id, "GammaLaw",     VAR_DBL, &(par->gammalaw), 0);
    io_write_attr_hdf5(id, "IO",     VAR_INT, &(par->io), 0);
    io_write_attr_hdf5(id, "NumCheckpoints",     VAR_INT, &(par->nChkpt), 0);
    io_write_attr_hdf5(id, "Init",         VAR_INT, &(par->init), 0);
    io_write_attr_hdf5(id, "InitPar1",     VAR_DBL, &(par->initPar1), 0);
    io_write_attr_hdf5(id, "InitPar2",     VAR_DBL, &(par->initPar2), 0);
    io_write_attr_hdf5(id, "InitPar3",     VAR_DBL, &(par->initPar3), 0);
    io_write_attr_hdf5(id, "InitPar4",     VAR_DBL, &(par->initPar4), 0);
    io_write_attr_hdf5(id, "InitPar5",     VAR_DBL, &(par->initPar5), 0);
    io_write_attr_hdf5(id, "InitPar6",     VAR_DBL, &(par->initPar6), 0);
    io_write_attr_hdf5(id, "InitPar7",     VAR_DBL, &(par->initPar7), 0);
    io_write_attr_hdf5(id, "InitPar8",     VAR_DBL, &(par->initPar8), 0);
}

void io_print_grid_hdf5(struct grid *g, struct parList *par, hid_t id)
{
    hsize_t dims1[1] = {1};
    hsize_t dims3[3] = {1,1,1};

    int nx1 = g->nx1;
    int nx2 = g->nx2;
    int nq = g->nq;

    dims1[0] = 1;
    io_write_dset_hdf5(id, "t", VAR_DBL, &(g->t), 1, dims1);
    dims1[0] = nx1+1;
    io_write_dset_hdf5(id, "x1", VAR_DBL, g->x1, 1, dims1);
    dims1[0] = nx2+1;
    io_write_dset_hdf5(id, "x2", VAR_DBL, g->x2, 1, dims1);

    dims3[0] = nx1;
    dims3[1] = nx2;
    dims3[2] = nq;
    io_write_dset_hdf5(id, "prim", VAR_DBL, g->prim, 3, dims3);
    io_write_dset_hdf5(id, "cons", VAR_DBL, g->cons, 3, dims3);
}

void io_write_attr_hdf5(hid_t id, char name[], int type, void *val, int len)
{
    /* 
     * Writes attribute to given hdf5 id. 'len' ignored unless type == VAR_STR.
     */

    hid_t space, dattype, attr;
    hsize_t dims[1] = {1};
    herr_t status;
    

    if(type == VAR_DBL)
        dattype = H5Tcopy(H5T_NATIVE_DOUBLE);
    else if(type == VAR_INT)
        dattype = H5Tcopy(H5T_NATIVE_INT);
    else if(type == VAR_LON)
        dattype = H5Tcopy(H5T_NATIVE_LONG);
    else if(type == VAR_STR)
    {
        dattype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(dattype, len);
    }
    else
        return;

    space = H5Screate_simple(1, dims, NULL);
    attr = H5Acreate(id, name, dattype, space, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr, dattype, val);

    status = H5Aclose(attr);
    status = H5Sclose(space);
    status = H5Tclose(dattype);
}

void io_write_dset_hdf5(hid_t id, char name[], int type, void *dat, int rank,
                        hsize_t dims[])
{
    /* 
     * Writes dataset to given hdf5 id.
     */

    hid_t space, dattype, dset;
    herr_t status;
    

    if(type == VAR_DBL)
        dattype = H5Tcopy(H5T_NATIVE_DOUBLE);
    else if(type == VAR_INT)
        dattype = H5Tcopy(H5T_NATIVE_INT);
    else if(type == VAR_LON)
        dattype = H5Tcopy(H5T_NATIVE_LONG);
    else
        return;

    space = H5Screate_simple(rank, dims, NULL);
    dset = H5Dcreate(id, name, dattype, space, H5P_DEFAULT, H5P_DEFAULT, 
                        H5P_DEFAULT);
    status = H5Dwrite(dset, dattype, H5S_ALL, H5S_ALL, H5P_DEFAULT, dat);

    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Tclose(dattype);
}
