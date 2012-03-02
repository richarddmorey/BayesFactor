


double thetaLogLikeAR_trend(double theta, double *beta, double *X, double sig2, double *y, int N, int p, double alphaTheta, double betaTheta);

double sampThetaAR_trend(double theta, double *beta, double *X, double sig2, double *y, int N, int p, double alphaTheta, double betaTheta, double sdmet);

void gibbsTwoSampleAR_trend(double *y, int N, double *X, int p, double rscaleInt, double rscaleSlp, double alphaTheta, double betaTheta, int iterations, double sdmet, double *chains, int progress, SEXP pBar, SEXP rho);

SEXP RgibbsTwoSampleAR_trend(SEXP yR, SEXP NR, SEXP XR, SEXP pR, SEXP rscaleIntR, SEXP rscaleSlpR, SEXP alphaThetaR, SEXP betaThetaR, SEXP iterationsR, SEXP sdmet, SEXP progressR, SEXP pBar, SEXP rho);