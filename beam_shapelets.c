// Original code by Kevin Huffenberger (2008)
#include <beam_shapelets.h>

void beam_shapelets_vec(double *val,double *rxvec, double *ryvec, int Npt, double beamx, double beamy, double gam, double fwhm_mean, double axisrat, double *shap_coef, int nmax);

double beam_shapelets(double rx, double ry, double beamx, double beamy, double gam, double fwhm_mean, double axisrat, double *shap_coef, int nmax){
  // Compute a single point in a beam, composed of a sum of beam shapelets

  double val;
  int Npt = 1;

  beam_shapelets_vec(&val,&rx, &ry, Npt, beamx, beamy, gam, fwhm_mean, axisrat, shap_coef,  nmax);

  return(val);
}



void beam_shapelets_vec(double *val,double *rxvec, double *ryvec, int Npt, double beamx, double beamy, double gam, double fwhm_mean, double axisrat, double *shap_coef, int nmax){
  //  Compute the beam as a sum of shapelets


  double sgam = sin(gam/180*M_PI);
  double cgam = cos(gam/180*M_PI);
  
  double fwhm_x = fwhm_mean/sqrt(axisrat);

  double sig_x = fwhm_x /sqrt(8*log(2));
  double sig_y = sig_x*axisrat;
 
  double u,v,x,y;

  int ipt;
  double rx,ry;
  double C1[nmax+1],C2[nmax+1];
  double Gx[nmax+1], Gy[nmax+1];
  int n1,n2,idx;
  
  
  shapelet_precompute_normHermite_coef(C1, C2, nmax);


  // loop over poiints
  for (ipt=0;ipt<Npt;ipt++) {
    rx = rxvec[ipt];
    ry = ryvec[ipt];
    val[ipt] = 0.0;
    
    // Scale and rotate into standard reference frame
    u = rx-beamx;
    v = ry-beamy;
    
    x = (cgam*u+sgam*v)/sig_x;
    y = (-sgam*u+cgam*v)/sig_y;
    
    
    // Compute and add in the shapelet components
    
    shapelet_normHermite(x,Gx,C1,C2,nmax);
    shapelet_normHermite(y,Gy,C1,C2,nmax);
    
    
    for (n1=0;n1<=nmax;n1++) {
      for (n2=0;n2+n1<=nmax;n2++) {
	
	idx = shapelet_cart_idx(n1,n2);
	
	val[ipt] += shap_coef[idx] * shapelet_cart(x, y, Gx, Gy, n1, n2);
    
	//	printf("%d %e %e\n",idx,shap_coef[idx],val[ipt]);
      }
    }
  }



}

void beam_shapelets_set(double *val,double *rxvec, double *ryvec, int Npt, double beamx, double beamy, double gam, double fwhm_mean, double axisrat, double *shap_coef, int nmax){
  // Compure the whole set of beam shapelets up to nmax, keeping them separate
  //
  // val needs size Nshapelet(nmax)*Npt

  double sgam = sin(gam/180*M_PI);
  double cgam = cos(gam/180*M_PI);
  
  double fwhm_x = fwhm_mean/sqrt(axisrat);

  double sig_x = fwhm_x /sqrt(8*log(2));
  double sig_y = sig_x*axisrat;
 
  double u,v,x,y;

  int ipt;
  double rx,ry;
  double C1[nmax+1],C2[nmax+1];
  double Gx[nmax+1], Gy[nmax+1];
  int n1,n2,idx;
  
  // printf("sig_x = %e\n",sig_x);
  
  shapelet_precompute_normHermite_coef(C1, C2, nmax);


  // loop over poiints
  for (ipt=0;ipt<Npt;ipt++) {
    rx = rxvec[ipt];
    ry = ryvec[ipt];
    
    // Scale and rotate into standard reference frame
    u = rx-beamx;
    v = ry-beamy;
    
    x = (cgam*u+sgam*v)/sig_x;
    y = (-sgam*u+cgam*v)/sig_y;
    
    
    // Compute and add in the shapelet components
    
    shapelet_normHermite(x,Gx,C1,C2,nmax);
    shapelet_normHermite(y,Gy,C1,C2,nmax);
    
    
    for (n1=0;n1<=nmax;n1++) {
      for (n2=0;n2+n1<=nmax;n2++) {
	
	idx = shapelet_cart_idx(n1,n2);
	
	val[idx*Npt+ipt] = shap_coef[idx] * shapelet_cart(x, y, Gx, Gy, n1, n2);
    
	//	printf("%d %e %e\n",idx,shap_coef[idx],val[ipt]);
      }
    }
  }



}



void beam_shapelets_evalset(double *eval, double *shapelets, int Npt, double *shap_coef, int nmax){
  // Evaluate a beam from precomputed shapelets on a given point set.

  int Nshapelet = shapelet_Ncoef(nmax); 

  int ish, ipt;
  
  for (ipt=0;ipt<Npt;ipt++) eval[ipt] = 0;
  
  for (ish=0;ish<Nshapelet;ish++) {
    for (ipt=0;ipt<Npt;ipt++) {
      eval[ipt] += shapelets[ish*Npt+ipt] * shap_coef[ish];
    }
  }

}


void beam_shapelets_one(double *val,double *rxvec, double *ryvec, int Npt, double beamx, double beamy, double gam, double fwhm_mean, double axisrat, double shap_coef, int n1, int n2, int nmax){
  //  Compute a map of a single shapelet


  double sgam = sin(gam/180*M_PI);
  double cgam = cos(gam/180*M_PI);
  
  double fwhm_x = fwhm_mean/sqrt(axisrat);

  double sig_x = fwhm_x /sqrt(8*log(2));
  double sig_y = sig_x*axisrat;
 
  double u,v,x,y;

  int ipt;
  double rx,ry;
  double C1[nmax+1],C2[nmax+1];
  double Gx[nmax+1], Gy[nmax+1];
  
  
  shapelet_precompute_normHermite_coef(C1, C2, nmax);


  // loop over poiints
  for (ipt=0;ipt<Npt;ipt++) {
    rx = rxvec[ipt];
    ry = ryvec[ipt];
    
    // Scale and rotate into standard reference frame
    u = rx-beamx;
    v = ry-beamy;
    
    x = (cgam*u+sgam*v)/sig_x;
    y = (-sgam*u+cgam*v)/sig_y;
    
    
    // Compute and add in the shapelet components
    
    shapelet_normHermite(x,Gx,C1,C2,nmax);
    shapelet_normHermite(y,Gy,C1,C2,nmax);
    
    
    val[ipt] = shap_coef * shapelet_cart(x, y, Gx, Gy, n1, n2);
    
  }


}
