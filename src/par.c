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
        return 0;

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

void print_pars(struct parList *theParList)
{
    printf("===Input Parameters===\n");
    printf("Hydro: %d\n", theParList->hydro);
    printf("EOS: %d\n", theParList->eos);
    printf("Cool: %d\n", theParList->cool);
    printf("Recon: %d\n", theParList->recon);
    printf("Riemann: %d\n", theParList->riemann);
    printf("Timestep: %d\n", theParList->step);
    printf("Frame: %d\n", theParList->frame);
    printf("BCX1Inner: %d\n", theParList->bc1Inner);
    printf("BCX1Outer: %d\n", theParList->bc1Outer);
    printf("BCX2Inner: %d\n", theParList->bc2Inner);
    printf("BCX2Outer: %d\n", theParList->bc2Outer);
    printf("Nx1: %d\n", theParList->nx1);
    printf("Nx2: %d\n", theParList->nx2);
    printf("Nghost: %d\n", theParList->nghost);
    printf("Ncons: %d\n", theParList->nc);
    printf("Npass: %d\n", theParList->np);
    printf("Tmin: %g\n", theParList->tmin);
    printf("Tmax: %g\n", theParList->tmax);
    printf("X1min: %g\n", theParList->x1min);
    printf("X1max: %g\n", theParList->x1max);
    printf("X2min: %g\n", theParList->x2min);
    printf("X2max: %g\n", theParList->x2max);
    printf("PLM: %g\n", theParList->plm);
    printf("CFL: %g\n", theParList->cfl);
    printf("GammaLaw: %g\n", theParList->gammalaw);
    printf("M: %g\n", theParList->M);
    printf("NumCheckpoints: %d\n", theParList->nChkpt);
    printf("Init: %d\n", theParList->init);
    printf("InitPar1: %g\n", theParList->initPar1);
    printf("InitPar2: %g\n", theParList->initPar2);
    printf("InitPar3: %g\n", theParList->initPar3);
    printf("InitPar4: %g\n", theParList->initPar4);
    printf("InitPar5: %g\n", theParList->initPar5);
    printf("InitPar6: %g\n", theParList->initPar6);
    printf("InitPar7: %g\n", theParList->initPar7);
    printf("InitPar8: %g\n", theParList->initPar8);
    printf("\n");
}
