#ifndef HYDRO2D_GRID
#define HYDRO2D_GRID

struct parList;

struct grid
{
    int nx1;
    int nx2;
    int nx1_int;
    int nx2_int;
    int ng11;
    int ng12;
    int ng21;
    int ng22;
    int d1;
    int d2;
    int nc;
    int nq;
    double x1min;
    double x1max;
    double x2min;
    double x2max;
    double t;
    double *x1;
    double *x2;
    double *prim;
    double *cons;
    double *cons_rk;
    double *prim_grad;
    double *prim_grad2;
    double PLM;
};

const static struct grid GRID_DEFAULT = {
    .nx1 = 0,
    .nx2 = 0,
    .nx1_int = 0,
    .nx2_int = 0,
    .ng11 = 0,
    .ng12 = 0,
    .ng21 = 0,
    .ng22 = 0,
    .d1 = 0,
    .d2 = 0,
    .nc = 0,
    .nq = 0,
    .x1min = 0.0,
    .x1max = 1.0,
    .x2min = 0.0,
    .x2max = 1.0,
    .t = 0.0,
    .x1 = NULL,
    .x2 = NULL,
    .prim = NULL,
    .prim_grad = NULL,
    .prim_grad2 = NULL,
    .cons = NULL,
    .cons_rk = NULL,
    .PLM = 1.5
};

void (*reconstruction)(struct grid *, int i, int j, int dir,
                        double gradlim[], double gradraw[], struct parList *par);

int set_reconstruction(struct parList *pars);
void make_grid(struct grid *g, struct parList *pars);
void free_grid(struct grid *g);

void calc_grad(struct grid *g, struct parList *par);
void interpolate_constant(struct grid *g, int i, int j, int dir,
                        double gradlim[], double gradraw[], struct parList *par);
void interpolate_plm(struct grid *g, int i, int j, int dir,
                        double gradlim[], double gradraw[], struct parList *par);

void copy_to_rk(struct grid *g);
void update_cons(struct grid *g, double fac1, double fac2);
void update_cons_rk(struct grid *g, double fac1, double fac2);
void update_x(struct grid *g, double fac1, double fac2);
void update_x_rk(struct grid *g, double fac1, double fac2);

#endif
