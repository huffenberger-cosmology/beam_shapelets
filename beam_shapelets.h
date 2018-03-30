#include <stdio.h>
#include <stdlib.h>
#include <math.h>




#include <shapelet_cartesian.h>


void beam_shapelets_set(double *val,double *rxvec, double *ryvec, int Npt, double beamx, double beamy, double gam, double fwhm_mean, double axisrat, double *shap_coef, int nmax);
double beam_shapelets(double rx, double ry, double beamx, double beamy, double gam, double fwhm_mean, double axisrat, double *shap_coef, int nmax);

void beam_shapelets_evalset(double *eval, double *shapelets, int Npt, double *shap_coef, int nmax);
