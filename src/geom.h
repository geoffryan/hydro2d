#ifndef HYDRO2D_GEOM
#define HYDRO2D_GEOM

struct parList;

void (*geom_CM)(double xm[], double xp[], double xc[]);
double (*geom_dA)(double xm[], double xp[], int dir);
double (*geom_dV)(double xm[], double xp[]);
void (*geom_gam)(double x[], double gam[3][3]);
void (*geom_igam)(double x[], double igam[3][3]);
void (*geom_dgam)(double x[], double dgam[2][3][3]);

int set_geometry(struct parList *par);

void geom_CM_flat_xyz(double x1[], double x2[], double xc[]);
double geom_dA_flat_xyz(double x1[], double x2[], int dir);
double geom_dV_flat_xyz(double x1[], double x2[]);
void geom_gam_flat_xyz(double x[], double gam[3][3]);
void geom_igam_flat_xyz(double x[], double gam[3][3]);
void geom_dgam_flat_xyz(double x[], double gam[2][3][3]);

#endif
