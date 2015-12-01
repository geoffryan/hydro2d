#ifndef HYDRO2D_PAR
#define HYDRO2D_PAR

enum{VAR_DBL, VAR_INT, VAR_LON, VAR_STR};

struct parList
{
    int hydro;
    int geom;
    int eos;
    int cool;
    int recon;
    int riemann;
    int step;
    int bc1Inner;
    int bc1Outer;
    int bc2Inner;
    int bc2Outer;
    int nx1;
    int nx2;
    int nghost;
    int nc;
    int np;
    double tmin;
    double tmax;
    double x1min;
    double x1max;
    double x2min;
    double x2max;
    double plm;
    double cfl;
    double gammalaw;

    int io;
    int nChkpt;

    int init;
    double initPar1;
    double initPar2;
    double initPar3;
    double initPar4;
    double initPar5;
    double initPar6;
    double initPar7;
    double initPar8;
};

const static struct parList PAR_DEFAULT = {
    .hydro = 0,
    .geom = 0,
    .eos = 0,
    .cool = 0,
    .recon = 0,
    .riemann = 0,
    .step = 0,
    .bc1Inner = 0,
    .bc1Outer = 0,
    .bc2Inner = 0,
    .bc2Outer = 0,
    .nx1 = 1,
    .nx2 = 1,
    .nghost = 1,
    .nc = 3,
    .np = 0,
    .tmin = 0.0,
    .tmax = 1.0,
    .x1min = 0.0,
    .x1max = 1.0,
    .x2min = 0.0,
    .x2max = 1.0,
    .plm = 1.5,
    .cfl = 0.5,
    .gammalaw = 1.4,

    .io = 0,
    .nChkpt = 0,

    .init = 0,
    .initPar1 = 0.0,
    .initPar2 = 0.0,
    .initPar3 = 0.0,
    .initPar4 = 0.0,
    .initPar5 = 0.0,
    .initPar6 = 0.0,
    .initPar7 = 0.0,
    .initPar8 = 0.0
};

int readvar(char filename[], char key[], int vtype, void *ptr);
void read_pars(struct parList *theParList, char filename[]);
void print_pars(struct parList *theParList, char filename[]);

#endif
