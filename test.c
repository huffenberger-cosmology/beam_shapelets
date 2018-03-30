#include <stdio.h>
#include <stdlib.h>

#include <decompose.h>
#include <shapelet_cartesian.h>

int main() {


  int Nx = 200;
  int Ny = 200;
  int Npt = Nx*Ny;

  int nmax = 5; // max shapelet order
  int Nshapelets = shapelet_Ncoef(nmax); // How many coefs


  double *shap_coef_in = calloc(Nshapelets,sizeof(double));  // input coefficients
  double *shap_coef_out = calloc(Nshapelets,sizeof(double));  // output coefficients
  double *shap_coef_unb = calloc(Nshapelets,sizeof(double));  // unbiased coefficients due to sampling correction

  double *shapelets = calloc(Nshapelets*Npt,sizeof(double));  // the shapelet functions
   
  double *overlap = calloc(Nshapelets*Nshapelets,sizeof(double)); // overlap integral and its inverse for debiasing
  double *overlapinv = calloc(Nshapelets*Nshapelets,sizeof(double));

  // coordinates for the "mother" gaussian beam function
  double beamx = 0.0;
  double beamy = 0.0;
  double gam = 30.0 * 180/M_PI; // orientation 
  double axisrat = 1.;  // ratio of semimajor axes
  double fwhm_mean = 1.0*sqrt(8*log(2)); // size in whatever units are being used (geometric mean of axes)

  double *shap_weight = decompose_shapelet_weight(Nshapelets, fwhm_mean);

  // Generate a set of points to evaluate the shapelets on
  // (here it's a grid but it's not required to be)
  int i,j,p;
  double minx = -5.;
  double miny = -5.;
  double maxx = 5.;
  double maxy = 5.;

  double Area = (maxx - minx)*(maxy-miny);
  
  double dx = (maxx-minx)/Nx;
  double dy = (maxy-miny)/Ny;

  double *ptx = calloc(Npt,sizeof(double));
  double *pty = calloc(Npt,sizeof(double));

  for (j=0;j<Ny;j++) {
    for (i=0;i<Nx;i++) {
      p = j*Nx + i;

      ptx[p] = minx + i*dx;
      pty[p] = miny + j*dy;
       
    }
  }
  
  // Evaluate the functions
  beam_shapelets_set(shapelets, ptx, pty, Npt, beamx, beamy, gam, fwhm_mean, axisrat, shap_weight, nmax);

  
  // Synthesize a test beam from coefficients
  double *testbeam = calloc(Npt,sizeof(double));

  shap_coef_in[shapelet_cart_idx(0,0)] = 1.0;
  shap_coef_in[shapelet_cart_idx(0,4)] = 0.25;
  shap_coef_in[shapelet_cart_idx(2,0)] = 0.1;
  
  beam_shapelets_evalset(testbeam, shapelets, Npt, shap_coef_in, nmax);

  // Decompose the test into coefficients
  double weight = Area/Npt;
  int ish;

  decompose_shapelets_set(Npt, testbeam, weight,  shapelets, Nshapelets, shap_coef_out); 
  for (ish=0;ish<Nshapelets;ish++) {
    shap_coef_out[ish] /= shap_weight[ish];  // correct for the mother beam size
  }

  // Compute the overlap integral and inverse
  decompose_overlap_integral(Npt, shapelets, Nshapelets, Area, shap_weight[0], overlap);
  fprintf(stdout,"condition = %e\n",decompose_condition(Nshapelets, overlap));
  decompose_overlap_invert(Nshapelets, overlap, overlapinv);

  // debias the coefficients
  decompose_unbias(Nshapelets, overlapinv, shap_coef_out, shap_coef_unb);
 
  // Check the output
  for (ish=0;ish<Nshapelets;ish++) {
    printf("% e % e % e | % e\n",
	   
	   shap_coef_in[ish],
	   shap_coef_out[ish],
	   shap_coef_unb[ish],
	   shap_weight[ish]);
      }


  // dump the test beam to a file
  FILE *fp = fopen("testbeam.dat","w");
  for (j=0;j<Ny;j++) {
    for (i=0;i<Nx;i++) {
      p = j*Nx + i;
      fprintf(fp,"% e ",testbeam[p]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  
  return(0);

}
