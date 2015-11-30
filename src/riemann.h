#ifndef HYDRO2D_RIEMANN
#define HYDRO2D_RIEMANN

void (*riemann_flux)(double primL[], double primR[], double F[], int nq,
                    double x[2], int dir, struct parList *pars);

int set_riemann_solver(struct parList *pars);

void lax_friedrichs_flux(double primL[], double primR[], double F[], int nq,
                            double x[2], int dir, struct parList *pars);
void hll_flux(double primL[], double primR[], double F[], int nq,
                 double x[2], int dir, struct parList *pars);
void hllc_flux(double primL[], double primR[], double F[], int nq,
                 double x[2], int dir, struct parList *pars);

#endif
