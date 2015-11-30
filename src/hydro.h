#ifndef HYDRO2D_HYDRO
#define HYDRO2D_HYDRO

void (*prim2cons)(double *prim, double *cons, double x[2], double dV,
                    struct parList *pars);
void (*cons2prim)(double *cons, double *prim, double x[2], double dV,
                    struct parList *pars);
void (*flux)(double *prim, double *F, double x[2], int dir, struct parList *pars);
void (*add_source)(double *prim, double *cons, double x[2], double dVdt, 
                    struct parList *pars);
void (*wave_speeds)(double *prim1, double *prim2, double *sL, double *sR,
                    double *sC, double x[2], int dir, struct parList *pars);
double (*mindt)(double *prim, double x[2], double dx[2],
                struct parList *pars);
void (*Ustar)(double *prim, double *Us, double sK, double sC, double x[2],
                    struct parList *pars);

int set_hydro(struct parList *pars);

void prim2cons_newt2(double *prim, double *cons, double x[2], double dV,
                            struct parList *pars);
void cons2prim_newt2(double *cons, double *prim, double x[2], double dV,
                            struct parList *pars);
void flux_newt2(double *prim, double *F, double x[2], int dir, 
                            struct parList *pars);
void add_source_newt2(double *prim, double *cons, double x[2], double dVdt, 
                            struct parList *pars);
void wave_speeds_newt2(double *prim1, double *prim2, double *sL, 
                            double *sR, double *sC, double x[2], int dir,
                            struct parList *pars);
double mindt_newt2(double *prim, double x[2], double dx[2], 
                            struct parList *pars);
void Ustar_newt2(double *prim, double *Us, double sK, double sC, 
                            double x[2], struct parList *pars);

#endif
