#ifndef HYDRO2D_INITIAL
#define HYDRO2D_INITIAL

struct grid;
struct parList;

void (*initial_value)(double *prim, double x[2], int nq, struct parList *par);

int set_initial(struct parList *par);
void initialize_grid(struct grid *g, struct parList *pars);

void initial_uniform(double *prim, double x[2], int nq, struct parList *par);

#endif
