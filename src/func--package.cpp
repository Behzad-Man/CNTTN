# include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

/// Check if simpler subalgorithm is appropriate.
inline bool CheckSimple(const double low, ///< lower bound of distribution
const double high ///< upper bound of distribution
) {
// Init Values Used in Inequality of Interest
double val1 = (2 * sqrt(exp(1))) / (low + sqrt(pow(low, 2) + 4));
double val2 = exp((pow(low, 2) - low * sqrt(pow(low, 2) + 4)) / (4)) ;
//

// Test if Simple is Preferred
if (high > low + val1 * val2) {
return true ;
} else {
return false ;
}
}

/// Draw using algorithm 1.

///
/// Naive Accept-Reject algorithm.
///
inline double UseAlgn1(const double low, ///< lower bound of distribution
const double high ///< upper bound of distribution
) {
// Init Valid Flag
int valid = 0 ;
//

// Init Draw Storage
double z = 0.0 ;
//

// Loop Until Valid Draw
while (valid == 0) {
z = Rf_rnorm(0.0, 1.0) ;

if (z <= high && z >= low) {
valid = 1 ;
}
}
//

// Returns
return z ;
//
}

/// Draw using algorithm 2.

///
///  Accept-Reject Algorithm
///

inline double UseAlgn2(const double low ///< lower bound of distribution
) {
// Init Values
const double alphastar = (low +
sqrt(pow(low, 2) + 4.0)
) / (2.0) ;
const double alpha = alphastar ;
double e = 0 ;
double z = 0 ;
double rho = 0 ;
double u = 0 ;
//

// Init Valid Flag
int valid = 0 ;
//

// Loop Until Valid Draw
while (valid == 0) {
e = Rf_rexp(1.0) ;
z = low + e / alpha ;

rho = exp(-pow(alpha - z, 2) / 2) ;
u = Rf_runif(0, 1) ;
if (u <= rho) {
// Keep Successes
valid = 1 ;
}
}
//

// Returns
return z ;
//
}

/// Draw using algorithm 3.

///
/// Accept-Reject Algorithm
///

inline double UseAlgn3(const double low, ///< lower bound of distribution
const double high ///< upper bound of distribution
) {
// Init Valid Flag
int valid = 0 ;
//

// Declare Qtys
double rho = 0 ;
double z = 0 ;
double u = 0 ;
//

// Loop Until Valid Draw
while (valid == 0) {
z = Rf_runif(low, high) ;
if (0 < low) {
rho = exp((pow(low, 2) - pow(z, 2)) / 2) ;
} else if (high < 0) {
rho = exp((pow(high, 2) - pow(z, 2)) / 2) ;
} else if (0 < high && low < 0) {
rho = exp(- pow(z, 2) / 2) ;
}

u = Rf_runif(0, 1) ;
if (u <= rho) {
valid = 1 ;
}
}
//

// Returns
return z ;
//
}


double rrtn1(const double mean,
const double sd,
const double low,
const double high
) {
// Namespace
using namespace Rcpp ;
//
// Init Useful Values
double draw = 0;
int type = 0 ;
int valid = 0 ; // used only when switching to a simplified version
// of Alg 2 within Type 4 instead of the less
// efficient Alg 3
//

// Set Current Distributional Parameters
const double c_mean = mean ;
double c_sd = sd ;
const double c_low = low ;
const double c_high = high ;
double c_stdlow = (c_low - c_mean) / c_sd ;
double c_stdhigh = (c_high - c_mean) / c_sd ; // bounds are standardized
//

// Map Conceptual Cases to Algorithm Cases
// Case 1 (Simple Deterministic AR)
if (0 <= c_stdhigh &&
0 >= c_stdlow
) {
type = 1 ;
}

// Case 2 (Robert 2009 AR)
// mu < low, high = Inf
if (0 < c_stdlow &&
c_stdhigh == INFINITY
) {
type = 2 ;
}

// Case 3 (Robert 2009 AR)
// high < mu, low = -Inf
if (0 > c_stdhigh &&
c_stdlow == -INFINITY
) {
type = 3 ;
}

// Case 4 (Robert 2009 AR)
if ((0 > c_stdhigh || 0 < c_stdlow) &&
!(c_stdhigh == INFINITY || c_stdlow == -INFINITY)
) {
type = 4 ;
}

////////////
// Type 1 //
////////////
if (type == 1) {
draw = UseAlgn1(c_stdlow, c_stdhigh) ;
}

////////////
// Type 3 //
////////////
if (type == 3) {
c_stdlow = -1 * c_stdhigh ;
c_stdhigh = INFINITY ;
c_sd = -1 * c_sd ; // hack to get two negative signs to cancel out

// Use Algorithm #2 Post-Adjustments
type = 2 ;
}

////////////
// Type 2 //
////////////
if (type == 2) {
draw = UseAlgn2(c_stdlow) ;
}

////////////
// Type 4 //
////////////
if (type == 4) {
if (CheckSimple(c_stdlow, c_stdhigh)) {
while (valid == 0) {
draw = UseAlgn2(c_stdlow) ;
// use the simple
// algorithm if it is more
// efficient
if (draw <= c_stdhigh) {
valid = 1 ;
}
}
} else {
draw = UseAlgn3(c_stdlow, c_stdhigh) ; // use the complex
// algorithm if the simple
// is less efficient
}
}



// Returns
return  c_mean + c_sd * draw ;
//
}
// [[Rcpp::export]]
SEXP rcpptno(SEXP sigma, SEXP a, SEXP mu, SEXP sig, SEXP ssig, SEXP low, SEXP high)
{
arma::mat H = Rcpp::as<arma::mat>(sigma);
arma::rowvec x = Rcpp::as<arma::rowvec>(a);
arma::rowvec u = Rcpp::as<arma::rowvec>(mu);
arma::rowvec s = Rcpp::as<arma::rowvec>(sig);
arma::rowvec ss = Rcpp::as<arma::rowvec>(ssig);
arma::rowvec alow = Rcpp::as<arma::rowvec>(low);
arma::rowvec bhigh = Rcpp::as<arma::rowvec>(high);
int m = x.size();
arma::rowvec ysim = x;
arma::rowvec V(m);
GetRNGstate();
for (int row=0; row<m; row++) {
V = ysim - u;
V[row] = 0;
ysim[row] = rrtn1(u[row] - s[row]*arma::as_scalar(H.row(row)*trans(V)),ss[row],
alow[row],bhigh[row]);
}
PutRNGstate();
return Rcpp::wrap(ysim);
}
// [[Rcpp::export]]
SEXP rcpptnm(SEXP sigma, SEXP a, SEXP mu, SEXP sig, SEXP ssig, SEXP low, SEXP high)
{
arma::mat H = Rcpp::as<arma::mat>(sigma);
arma::mat x = Rcpp::as<arma::mat>(a);
arma::mat u = Rcpp::as<arma::mat>(mu);
arma::rowvec s = Rcpp::as<arma::rowvec>(sig);
arma::rowvec ss = Rcpp::as<arma::rowvec>(ssig);
arma::mat alow = Rcpp::as<arma::mat>(low);
arma::mat bhigh = Rcpp::as<arma::mat>(high);
int n = x.n_cols;
int m = ss.size();
arma::mat ysim = x;
arma::colvec V(m);
GetRNGstate();
for (int col=0; col<n; col++) {
for (int row=0; row<m; row++) {
V = ysim.col(col) - u.col(col);
V[row] = 0;
ysim(row,col) = rrtn1(u(row,col) - s[row]*arma::as_scalar(H.row(row)*V),ss[row],alow(row,col),bhigh(row,col));
}
}
PutRNGstate();
return Rcpp::wrap(ysim);
}
// [[Rcpp::export]]
SEXP rcppcomdel(SEXP W, SEXP U, SEXP H)
{
arma::mat Wc = Rcpp::as<arma::mat>(W);
arma::mat Uc = Rcpp::as<arma::mat>(U);
arma::mat Hc = Rcpp::as<arma::mat>(H);
int n = Wc.n_cols;
arma::mat D = diagmat(Uc.col(0));
arma::mat Z = D*Hc;
arma::colvec V = Z*Wc.col(0);
arma::mat M = Z*D;
for (int col=1; col<n; col++) {
D = diagmat(Uc.col(col));
Z = D*Hc;
V = Z*Wc.col(col) + V;
M = Z*D + M;
}
return Rcpp::List::create(Rcpp::Named("termd1") = V,Rcpp::Named("termd2") = M);
}

// [[Rcpp::export]]
SEXP rcppng(SEXP sigma, SEXP sigw, SEXP isigt, SEXP invL, SEXP muW, SEXP rngN)
{
  arma::mat H = Rcpp::as<arma::mat>(sigma);
  arma::mat x = Rcpp::as<arma::mat>(invL);
  arma::mat mu = Rcpp::as<arma::mat>(muW);
  arma::mat nsim = Rcpp::as<arma::mat>(rngN);
  double sgww = Rcpp::as<double>(sigw);
  double ittau = Rcpp::as<double>(isigt);
  int n = x.n_cols;
  int m = x.n_rows;
  arma::mat V(m,n);
    for (int col=0; col<n; col++) {
    arma::mat KW = inv_sympd (H + sgww*ittau*diagmat(x.col(col)));
    V.col(col) = KW*mu.col(col)+chol(KW,"lower")*nsim.col(col);
  }
  return Rcpp::wrap(V);
}
