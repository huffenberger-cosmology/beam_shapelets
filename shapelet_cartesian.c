// Original code by Kevin Huffenberger (2008)
#include <math.h>
#include <shapelet_cartesian.h>


#define SQRT_PI 1.7724538509055158819194275565678

int shapelet_cart_idx(int n1,int n2) {

  return( (n1+n2)*((n1+n2)+1)/2 + n2);

}

void shapelet_cart_2idx(int idx, int *n1, int *n2) {
  int na=0, nb=0;
  
  
  na = (sqrt(8*idx + 1) - 1)/2;
  nb = idx - shapelet_cart_idx(na,0);


  *n1 = na-nb;
  *n2 = nb;

}

int shapelet_Ncoef(int nmax) {
  return( shapelet_cart_idx(nmax+1,0)  );
}

double shapelet_cart(double x, double y, double *Gx, double *Gy, int n1, int n2) {
  // Compute shapelet coefficient for beta = 1
  // Note, normalized to pi, not 1, so that the 0,0 component is 1 at peak.

  double phi;

  phi = Gx[n1] * Gy[n2] * exp(-(x*x + y*y)/2) / pow(2,(n1+n2)/2.);
  
  return(phi);
}



void shapelet_normHermite_old(double x, double *G, int nmax){
  //returns an array of G_n(x) = H_n(x)/sqrt(n!)
  // G is size nmax+1
  
  G[0] = 1;
  G[1] = 2*x;

  int n;
  for (n = 2;n <= nmax; n++) {
    G[n] = 2*x/sqrt(n)*G[n-1] - 2*sqrt((n-1)/(double) n)*G[n-2];
  }

}

void shapelet_precompute_normHermite_coef(double *C1, double *C2, int nmax){
  // precompute the coefficients for the hermite recursion:
  // arrays are for 
  //
  // C1 = 2/sqrt(n) 
  // C2 = 2 sqrt(n-1)/sqrt(n)
  //
  // and should be nmax+1 elements
  
  int n;
  for (n=0;n<=nmax;n++) {
    C1[n] = 2/sqrt(n);
    C2[n] = 2*sqrt((n-1)/(double) n);
  }
}  


void shapelet_normHermite(double x, double *G, double *C1, double *C2, int nmax) {
  //returns an array of G_n(x) = H_n(x)/sqrt(n!)
  // G is size nmax+1
  
  G[0] = 1;
  G[1] = 2*x;
  
  int n;
  for (n = 2;n <= nmax; n++) {
    G[n] = C1[n]*x*G[n-1] - C2[n]*G[n-2];
  }
  
  
}
  



