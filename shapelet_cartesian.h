int shapelet_cart_idx(int n1,int n2);

int shapelet_Ncoef(int nmax);
void shapelet_normHermite(double x, double *G, double *C1, double *C2, int nmax);

double shapelet_cart(double x, double y, double *Gx, double *Gy, int n1, int n2);
void shapelet_precompute_normHermite_coef(double *C1, double *C2, int nmax);
