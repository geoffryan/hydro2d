#ifndef HYDRO2D_BC
#define HYDRO2D_BC

struct grid;
struct parList;

void (*bc_1L)(struct grid *g, struct parList *par);
void (*bc_1R)(struct grid *g, struct parList *par);
void (*bc_2L)(struct grid *g, struct parList *par);
void (*bc_2R)(struct grid *g, struct parList *par);

int set_boundary(struct parList *par);

void bc_1L_none(struct grid *g, struct parList *par);
void bc_1R_none(struct grid *g, struct parList *par);
void bc_2L_none(struct grid *g, struct parList *par);
void bc_2R_none(struct grid *g, struct parList *par);

void bc_1L_fixed(struct grid *g, struct parList *par);
void bc_1R_fixed(struct grid *g, struct parList *par);
void bc_2L_fixed(struct grid *g, struct parList *par);
void bc_2R_fixed(struct grid *g, struct parList *par);

void bc_1L_outflow(struct grid *g, struct parList *par);
void bc_1R_outflow(struct grid *g, struct parList *par);
void bc_2L_outflow(struct grid *g, struct parList *par);
void bc_2R_outflow(struct grid *g, struct parList *par);

void bc_1L_periodic(struct grid *g, struct parList *par);
void bc_1R_periodic(struct grid *g, struct parList *par);
void bc_2L_periodic(struct grid *g, struct parList *par);
void bc_2R_periodic(struct grid *g, struct parList *par);

#endif
