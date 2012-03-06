


double thetaLogLikeAR_trend(double theta, double *beta, double *X, double sig2, double *y, int N, int p, double alphaTheta, double betaTheta);

double sampThetaAR_trend(double theta, double *beta, double *X, double sig2, double *y, int N, int p, double alphaTheta, double betaTheta, double sdmet);

void gibbsTwoSampleAR_trend(double *y, int N, double *X, int p, double rscaleInt, double rscaleSlp, double alphaTheta, double betaTheta,  double loInt, double upInt, double loSlp, double upSlp, int iterations, double sdmet, double *chains, double *postdens, int progress, SEXP pBar, SEXP rho);

SEXP RgibbsTwoSampleAR_trend(SEXP yR, SEXP NR, SEXP XR, SEXP pR, SEXP rscaleIntR, SEXP rscaleSlpR, SEXP alphaThetaR, SEXP betaThetaR, SEXP loIntR, SEXP upIntR, SEXP loSlpR, SEXP upSlpR, SEXP iterationsR, SEXP sdmet, SEXP progressR, SEXP pBar, SEXP rho);