#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <beam_shapelets.h>
#include <shapelet_cartesian.h>

#include <string.h>



double *decompose_shapelet_weight(int Nshapelets, double fwhm_mean);



void decompose_shapelets_set(int Npt, double *f, double weight,  double *shapelets, int Nshapelet, double *shap_coef_out);
void decompose_overlap_integral(int Npt, double *shapelets, int Nshapelet, double A, double shapelet_norm, double *overlap);
void decompose_overlap_invert(int Nshapelets, double *overlap, double *overlapinv);
double decompose_condition(int Nshapelet,  double *overlap);
void decompose_unbias(int Nshapelets, double *overlapinv, double *shap_coef_bias, double *shap_coef_unbias);
