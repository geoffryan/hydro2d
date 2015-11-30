#include <stdio.h>
#include <string.h>
#include "par.h"

int readvar(char filename[], char key[], int vtype, void *ptr)
{
    FILE *f = fopen(filename, "r");

    char line[256];
    char word[256];
    int found = 0;

    while(fgets(line,256,f) != NULL)
    {
        sscanf(line, "%s ", word);
        if(strcmp(word,key) == 0)
        {
            found = 1;
            break;
        }
    }
    fclose(f);
    if(!found)
    {
        printf("SETUP: %s parameter not found. Using default.\n", key);
        return 1;
    }

    char *sval = line + strlen(key) + strspn(line+strlen(key)," \t:=");

    if(vtype == VAR_DBL)
    {
        double val;
        sscanf(sval, "%lf", &val);
        *((double *)ptr) = val;
    }
    else if(vtype == VAR_INT)
    {
        int val;
        sscanf(sval, "%d", &val);
        *((int *)ptr) = val;
    }
    else if(vtype == VAR_LON)
    {
        long val;
        sscanf(sval, "%ld", &val);
        *((long *)ptr) = val;
    }
    else
    {
        strcpy((char *) ptr, sval);
    }

    return 0;
}

void read_pars(struct parList *theParList, char filename[])
{
    readvar(filename, "Hydro", VAR_INT, &(theParList->hydro));
    readvar(filename, "Geometry", VAR_INT, &(theParList->geom));
    readvar(filename, "EOS",   VAR_INT, &(theParList->eos));
    readvar(filename, "Cool",  VAR_INT, &(theParList->cool));
    readvar(filename, "Recon",  VAR_INT, &(theParList->recon));
    readvar(filename, "Riemann",  VAR_INT, &(theParList->riemann));
    readvar(filename, "Timestep",  VAR_INT, &(theParList->step));
    readvar(filename, "Frame",  VAR_INT, &(theParList->frame));
    readvar(filename, "BCX1Inner",  VAR_INT, &(theParList->bc1Inner));
    readvar(filename, "BCX1Outer",  VAR_INT, &(theParList->bc1Outer));
    readvar(filename, "BCX2Inner",  VAR_INT, &(theParList->bc2Inner));
    readvar(filename, "BCX2Outer",  VAR_INT, &(theParList->bc2Outer));
    readvar(filename, "Nx1",     VAR_INT, &(theParList->nx1));
    readvar(filename, "Nx2",     VAR_INT, &(theParList->nx2));
    readvar(filename, "Nghost",     VAR_INT, &(theParList->nghost));
    readvar(filename, "Ncons",     VAR_INT, &(theParList->nc));
    readvar(filename, "Npass",     VAR_INT, &(theParList->np));
    readvar(filename, "Tmin",     VAR_DBL, &(theParList->tmin));
    readvar(filename, "Tmax",     VAR_DBL, &(theParList->tmax));
    readvar(filename, "X1min",     VAR_DBL, &(theParList->x1min));
    readvar(filename, "X1max",     VAR_DBL, &(theParList->x1max));
    readvar(filename, "X2min",     VAR_DBL, &(theParList->x2min));
    readvar(filename, "X2max",     VAR_DBL, &(theParList->x2max));
    readvar(filename, "PLM",     VAR_DBL, &(theParList->plm));
    readvar(filename, "CFL",     VAR_DBL, &(theParList->cfl));
    readvar(filename, "GammaLaw",     VAR_DBL, &(theParList->gammalaw));
    readvar(filename, "M",     VAR_DBL, &(theParList->M));
    readvar(filename, "IO",     VAR_INT, &(theParList->io));
    readvar(filename, "NumCheckpoints",     VAR_INT, &(theParList->nChkpt));
    readvar(filename, "Init",         VAR_INT, &(theParList->init));
    readvar(filename, "InitPar1",     VAR_DBL, &(theParList->initPar1));
    readvar(filename, "InitPar2",     VAR_DBL, &(theParList->initPar2));
    readvar(filename, "InitPar3",     VAR_DBL, &(theParList->initPar3));
    readvar(filename, "InitPar4",     VAR_DBL, &(theParList->initPar4));
    readvar(filename, "InitPar5",     VAR_DBL, &(theParList->initPar5));
    readvar(filename, "InitPar6",     VAR_DBL, &(theParList->initPar6));
    readvar(filename, "InitPar7",     VAR_DBL, &(theParList->initPar7));
    readvar(filename, "InitPar8",     VAR_DBL, &(theParList->initPar8));
}

void print_pars(struct parList *theParList, char filename[])
{
    FILE *f;
    if(filename == NULL)
        f = stdout;
    else
        f = fopen(filename, "a");

    fprintf(f, "### Input Parameters ###\n");
    fprintf(f, "Hydro: %d\n", theParList->hydro);
    fprintf(f, "Geometry: %d\n", theParList->geom);
    fprintf(f, "EOS: %d\n", theParList->eos);
    fprintf(f, "Cool: %d\n", theParList->cool);
    fprintf(f, "Recon: %d\n", theParList->recon);
    fprintf(f, "Riemann: %d\n", theParList->riemann);
    fprintf(f, "Timestep: %d\n", theParList->step);
    fprintf(f, "Frame: %d\n", theParList->frame);
    fprintf(f, "BCX1Inner: %d\n", theParList->bc1Inner);
    fprintf(f, "BCX1Outer: %d\n", theParList->bc1Outer);
    fprintf(f, "BCX2Inner: %d\n", theParList->bc2Inner);
    fprintf(f, "BCX2Outer: %d\n", theParList->bc2Outer);
    fprintf(f, "Nx1: %d\n", theParList->nx1);
    fprintf(f, "Nx2: %d\n", theParList->nx2);
    fprintf(f, "Nghost: %d\n", theParList->nghost);
    fprintf(f, "Ncons: %d\n", theParList->nc);
    fprintf(f, "Npass: %d\n", theParList->np);
    fprintf(f, "Tmin: %g\n", theParList->tmin);
    fprintf(f, "Tmax: %g\n", theParList->tmax);
    fprintf(f, "X1min: %g\n", theParList->x1min);
    fprintf(f, "X1max: %g\n", theParList->x1max);
    fprintf(f, "X2min: %g\n", theParList->x2min);
    fprintf(f, "X2max: %g\n", theParList->x2max);
    fprintf(f, "PLM: %g\n", theParList->plm);
    fprintf(f, "CFL: %g\n", theParList->cfl);
    fprintf(f, "GammaLaw: %g\n", theParList->gammalaw);
    fprintf(f, "M: %g\n", theParList->M);
    fprintf(f, "IO: %d\n", theParList->io);
    fprintf(f, "NumCheckpoints: %d\n", theParList->nChkpt);
    fprintf(f, "Init: %d\n", theParList->init);
    fprintf(f, "InitPar1: %g\n", theParList->initPar1);
    fprintf(f, "InitPar2: %g\n", theParList->initPar2);
    fprintf(f, "InitPar3: %g\n", theParList->initPar3);
    fprintf(f, "InitPar4: %g\n", theParList->initPar4);
    fprintf(f, "InitPar5: %g\n", theParList->initPar5);
    fprintf(f, "InitPar6: %g\n", theParList->initPar6);
    fprintf(f, "InitPar7: %g\n", theParList->initPar7);
    fprintf(f, "InitPar8: %g\n", theParList->initPar8);

    if(filename != NULL)
        fclose(f);
}
