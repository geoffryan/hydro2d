#ifndef HYDRO2D_TIMESTEP
#define HYDRO2D_TIMESTEP

void (*timestep)(struct grid *, double, struct parList *);

int set_timestep(struct parList *pars);
double get_dt(struct grid *g, struct parList *pars);
void step_fe(struct grid *g, double dt, struct parList *pars);
void step_rk2_mp(struct grid *g, double dt, struct parList *pars);
void step_rk2_tvd(struct grid *g, double dt, struct parList *pars);
void step_rk3_tvd(struct grid *g, double dt, struct parList *pars);
void calc_cons(struct grid *g, struct parList *pars);
void calc_prim(struct grid *g, struct parList *pars);

void substep(struct grid *g, double rkfac1, double rkfac2, double dt, 
                struct parList *pars);

#endif
