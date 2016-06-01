#ifndef HYDRO2D_GEOM
#define HYDRO2D_GEOM

struct parList;

void (*geom_CM)(double xm[], double xp[], double xc[]);
void (*geom_CM2)(double xm[], double xp[], double xc[], int dir);
double (*geom_dA)(double xm[], double xp[], int dir);
double (*geom_dV)(double xm[], double xp[]);
double (*geom_J)(double x[]);
double (*geom_J2)(double x[], int dir);
void (*geom_gam)(double x[], double gam[3][3]);
void (*geom_igam)(double x[], double igam[3][3]);
void (*geom_dgam)(double x[], double dgam[2][3][3]);

int set_geometry(struct parList *par);

void geom_CM_cart_xyz(double x1[], double x2[], double xc[]);
void geom_CM2_cart_xyz(double x1[], double x2[], double xc[], int dir);
double geom_dA_cart_xyz(double x1[], double x2[], int dir);
double geom_dV_cart_xyz(double x1[], double x2[]);
double geom_J_cart_xyz(double x[]);
double geom_J2_cart_xyz(double x[], int dir);
void geom_gam_cart_xyz(double x[], double gam[3][3]);
void geom_igam_cart_xyz(double x[], double gam[3][3]);
void geom_dgam_cart_xyz(double x[], double gam[2][3][3]);

void geom_CM_cyl_rpz(double x1[], double x2[], double xc[]);
void geom_CM2_cyl_rpz(double x1[], double x2[], double xc[], int dir);
double geom_dA_cyl_rpz(double x1[], double x2[], int dir);
double geom_dV_cyl_rpz(double x1[], double x2[]);
double geom_J_cyl_rpz(double x[]);
double geom_J2_cyl_rpz(double x[], int dir);
void geom_gam_cyl_rpz(double x[], double gam[3][3]);
void geom_igam_cyl_rpz(double x[], double gam[3][3]);
void geom_dgam_cyl_rpz(double x[], double gam[2][3][3]);

void geom_CM_sph_tpr(double x1[], double x2[], double xc[]);
void geom_CM2_sph_tpr(double x1[], double x2[], double xc[], int dir);
double geom_dA_sph_tpr(double x1[], double x2[], int dir);
double geom_dV_sph_tpr(double x1[], double x2[]);
double geom_J_sph_tpr(double x[]);
double geom_J2_sph_tpr(double x[], int dir);
void geom_gam_sph_tpr(double x[], double gam[3][3]);
void geom_igam_sph_tpr(double x[], double gam[3][3]);
void geom_dgam_sph_tpr(double x[], double gam[2][3][3]);

void geom_CM_cyl_rpz_B(double x1[], double x2[], double xc[]);
void geom_CM2_cyl_rpz_B(double x1[], double x2[], double xc[], int dir);
double geom_dA_cyl_rpz_B(double x1[], double x2[], int dir);
double geom_dV_cyl_rpz_B(double x1[], double x2[]);
double geom_J_cyl_rpz_B(double x[]);
double geom_J2_cyl_rpz_B(double x[], int dir);
void geom_gam_cyl_rpz_B(double x[], double gam[3][3]);
void geom_igam_cyl_rpz_B(double x[], double gam[3][3]);
void geom_dgam_cyl_rpz_B(double x[], double gam[2][3][3]);

#endif
