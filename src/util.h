double beta_logpost(int &n, double *Y, double *eta, double phi);

double norm_logpost(int &n, double *Y, double *eta, double phi);

double tnorm_logpost(int &n, double *Y, double *eta, double phi, double a, double b);

void G(double *alpha, double *t, double *eta, int n);

double G(double *alpha, double t);

double mirror(double theta, double a, double b);

double sampleU(double mu, double tuning, double a, double b);

void zeros(double *a, int n);

double logit(double theta, double a, double b);

double logitInv(double z, double a, double b);

void printMtrx(double *m, int nRow, int nCol);

double dtnorm(double x, double mu, double sigma, double a, double b, int lg);

double yl(int k);

double rtexp(double a, double b);

double rtnorm(const double mu, const double sigma, double a, double b);
