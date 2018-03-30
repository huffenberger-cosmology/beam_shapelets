#include <decompose.h>


void decompose_shapelets(int Npt, double *rx, double *ry, double *f, double weight, double beamx, double beamy, double gam, double fwhm_mean, double axisrat, double *shap_coef_out, int nmax){
  
  
  int ipt;
  int n,np;
  int Nshapelets = shapelet_Ncoef(nmax); 

  double *shap_coef = calloc(Nshapelets,sizeof(double));;

  for (n=0;n<Nshapelets;n++) {
    shap_coef_out[n] = 0.0;
  }


  for (n=0;n<Nshapelets;n++) {
    
    for (np=0;np<Nshapelets;np++) {
      shap_coef[np] = 0;
    }
    shap_coef[n] = 1.0/ M_PI;



    for (ipt=0; ipt<Npt; ipt++) {
      
      shap_coef_out[n] +=  weight/pow(fwhm_mean/sqrt(8*log(2)),2)* f[ipt]* beam_shapelets(rx[ipt], ry[ipt], beamx, beamy, gam,  fwhm_mean, axisrat, shap_coef, nmax);
      //      shap_coef_out[n] += weight * pow(fwhm_mean,2) * f[ipt] * beam_shapelets(rx[ipt], ry[ipt], beamx, beamy, gam,  fwhm_mean, axisrat, shap_coef, nmax);

    }
  }

  free(shap_coef);
}


double *decompose_shapelet_weight(int Nshapelets, double fwhm_mean) {
  
  double *shap_weight = calloc(Nshapelets,sizeof(double));
  int ish;
  double sig = fwhm_mean/sqrt(8*log(2));
  double norm = 1.0/pow(sig,2)/M_PI;

  for (ish=0;ish<Nshapelets;ish++) {
    shap_weight[ish] = norm;
  }

  return(shap_weight);
}


void decompose_shapelets_set(int Npt, double *f, double weight,  double *shapelets, int Nshapelet, double *shap_coef_out){
  //  f has size Npt
  // shapelets has size Nshapelet * Npt
  // shap_coef_out has size Nshapelet


  int ipt;
  int n;
  
  for (n=0;n<Nshapelet;n++) {
    shap_coef_out[n] = 0.0;
  }
  
  for (n=0;n<Nshapelet;n++) {
    for (ipt=0; ipt<Npt; ipt++) {
      shap_coef_out[n] +=  f[ipt] * shapelets[n*Npt+ipt];
    }

    shap_coef_out[n] *= weight;
  }
  
}



void decompose_overlap_integral(int Npt, double *shapelets, int Nshapelet, double A, double shapelet_norm, double *overlap) {
  // shapelets has size Nshapelet * Npt
  // overlap has size Nshapelet * Nshapelet
  
  // shapelet_norm records how the shapelets are normalize...
  //           ... normally equal to  1.0/pow(fwhm_mean/sqrt(8*log(2)),2)/M_PI;
  // A is the bounding area of the pointing
  
  double norm = A/Npt/shapelet_norm;
  int ish,ish1;
  int ipt;




  for (ish1=0;ish1<Nshapelet;ish1++) {
    for (ish=0;ish<Nshapelet;ish++) {
      overlap[ish1*Nshapelet+ish] = 0;
      for (ipt=0;ipt<Npt;ipt++) {
	
	if (ish1>ish) overlap[ish1*Nshapelet+ish] = overlap[ish*Nshapelet+ish1];
	
	overlap[ish1*Nshapelet+ish] += norm * shapelets[ish*Npt+ipt]*shapelets[ish1*Npt+ipt];
	
      }
    }
  }
  

}

#ifndef CLAPACK
#define dpotrf dpotrf_
#define dpotri dpotri_
//      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
int dpotrf(char *UPLO, int *N, double *A, int *LDA, int *INFO);
//      SUBROUTINE DPOTRI( UPLO, N, A, LDA, INFO )
int dpotri(char *UPLO, int *N, double *A, int *LDA, int *INFO);
#else

int dpotrf(char UPLO, int N, double *A, int LDA, int *INFO);
int dpotri(char UPLO, int N, double *A, int LDA, int *INFO);
#endif

void decompose_overlap_invert(int Nshapelets, double *overlap, double *overlapinv) {

  int N =  Nshapelets;
  int N2 = N*N;
  double *U = calloc(N2,sizeof(double));
  int INFO = 0;
  char UPLO = 'U';

  memcpy(U,overlap,N2*sizeof(double));

#ifndef CLAPACK
  dpotrf(&UPLO, &N, U, &N, &INFO);
  dpotri(&UPLO, &N, U, &N, &INFO);
  //  printf("INFO = %d\n",INFO);
#else
  dpotrf(UPLO, N, U, N, &INFO);
  dpotri(UPLO, N, U, N, &INFO);
#endif  

  // Copy both halves into output
  int i,j;
  for (j=0;j<N;j++) {
    for (i=0;i<=j;i++) {
      overlapinv[j*N+i] = overlapinv[i*N+j] = U[j*N+i];
    }
  }

  free(U);
}

//      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
void dsyev_(char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO );

double decompose_condition(int Nshapelet,  double *overlap) {
  // Compute condition number
  char JOBZ = 'N';
  char UPLO = 'L';
  double *overlap_copy = calloc(Nshapelet*Nshapelet,sizeof(double));
  double *eigenvalues = calloc(Nshapelet,sizeof(double));
  double *WORK = calloc(4*Nshapelet,sizeof(double));
  int LWORK = 4*Nshapelet;
  int INFO = 0;
  
  memcpy(overlap_copy,overlap,Nshapelet*Nshapelet*sizeof(double));
  
  dsyev_(&JOBZ, &UPLO, &Nshapelet, overlap_copy, &Nshapelet, eigenvalues, WORK,
	&LWORK, &INFO );
  
  
  double condition = eigenvalues[Nshapelet-1]/eigenvalues[0];
  
  free(overlap_copy);
  free(eigenvalues);
  free(WORK);
  
  return(condition);
}




void decompose_unbias(int Nshapelets, double *overlapinv, double *shap_coef_bias, double *shap_coef_unbias) {

  int N = Nshapelets;
  int i,j;


  for (j=0;j<N;j++) {
    shap_coef_unbias[j] = 0.0;
    
    for (i=0;i<N;i++) {
      shap_coef_unbias[j] += overlapinv[j*N+i] * shap_coef_bias[i];
    }
  }
  
}

