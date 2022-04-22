/* Based on PointMetrics1D.c from Phillips 2015.  Biol Invasions 17:1946 

Estimates density of individuals, mean trait values, and trait variances at each individual's location
Assumes a Gaussian kernel as the neighbourhood size.
Requires an array of 1D locations (X[i]), trait values (H[i], D[i])
bandwidth currently fixed at bw=1 (easy enough to make variable)
 
 2D functions also follow at the end
*/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
//#include <Rdefines.h>
#include <math.h>
#include "Rinterface.h"
#include <stdbool.h>

#define Pi 3.141593
// #define x_max 50.0
// #define y_max 50.0


//-----------------------//
// Function declarations //
//-----------------------//
double norm1D (double x, double bw);
double norm2D (double z, double bw);
double ToroidalDistance (double x1, double y1, double x2, double y2, double x_min, double y_min, double x_max, double y_max);
double NonToroidalDistance (double x1, double y1, double x2, double y2);
SEXP metrics (SEXP R_X, SEXP R_Y, SEXP R_n, SEXP R_bw, SEXP R_torus, SEXP R_x_min, SEXP R_y_min, SEXP R_x_max, SEXP R_y_max);
SEXP mate (SEXP R_fX, SEXP R_fY, SEXP R_mX, SEXP R_mY, SEXP R_nF, SEXP R_nM, SEXP R_bw, SEXP R_torus, SEXP R_x_min, SEXP R_y_min, SEXP R_x_max, SEXP R_y_max);

//------------------------//
//  Function definitions  //
//------------------------//

// Probability density at x (location on a 1D line)
//		- given that x is normally distributed with mean 0 and sd=bw
double norm1D (double x, double bw){
	return exp(-pow(x,2)/(2*pow(bw,2)))/(bw*sqrt(2*Pi));
}

// Probability density at x and y (location on a 2D plane)
// Bivariate normal pdf = (e^(-z/(2*(1-covxy)))) / (2*pi*sdx*sdy*sqrt(1-covxy))
//		where: z = (((x-ux)^2)/sdx^2) - ((2*covxy*(x-ux)*(y-uy))/(sdx*sdy)) + (((y-uy)^2)/sdy^2)
//		   --> z = ((x-ux)^2 + (y-uy)^2) / sd^2
//		also:  c = sqrt(dx^2 + dy^2) is the distance between 2 points, dx=x2-x1, and dy=y2-y1
//		   --> c = sqrt(z * sd^2)
//		then:  z = (c^2) / (sd^2)
//		hence: pdf = (e^((-c^2)/(2*sd^2))) / (2*pi*sd^2);
//			since: covxy=0; mux=muy=0; sdx=sdy; z^2=x^2+y^2
double norm2D (double z, double bw){
	return exp(-pow(z,2)/(2*pow(bw,2)))/(2*Pi*pow(bw,2));
}

// Toroidal distance
//		- The landscape is an edgeless torus which wraps around on itself.
//		- This means that the maximum distance between 2 points along the x or y axes is half the width or height of the landscape, repectively.
//		- Specifially, if the distance between 2 points is greater than half of the landscape's,
//		- then its actual distance is the x_max - dx or y_max - dy.
double ToroidalDistance (double x1, double y1, double x2, double y2, double x_min, double y_min, double x_max, double y_max){
  double dx = fabs(x2 - x1);
  double dy = fabs(y2 - y1);
  if (dx > ((x_max+x_min)/2.0))
    dx = (x_max+x_min) - dx;
  if (dy > ((y_max+y_min)/2.0))
    dy = (y_max+y_min) - dy;
  return sqrt(dx*dx + dy*dy);
}

// Non-toroidal distance (flat bounded x_min to x_max and y_min to y_max)
double NonToroidalDistance (double x1, double y1, double x2, double y2){
  double dx = x2 - x1;
  double dy = y2 - y1;
  return sqrt(dx*dx + dy*dy);
}

// Calculate the probability density on a bivariate normal distribution given the
// pairwise toroidal or non-toroidal distances between individuals
// (i.e. distances on a toroidal surface or on a flat space, respectively).
// Specifically, calculate the density around eah individual ii
// by calculating the pairwise toroidal distance
// between individual ii and all the other individuals,
// and using it to calculate the density on a bivariate normal distribution.
// Hence, as the distance increases (z = sqrt(|x|^2+|y|^2)),
//			 the distance from the mean of the distribution also increases, and
//			 the density decreases.
// The density around individual ii is the sum of densities 
// across all pairwise distances with the other individuals.
// NOTE: for 1D simulations set y_min == y_max
SEXP metrics (SEXP R_X, SEXP R_Y, SEXP R_n, SEXP R_bw, SEXP R_torus, SEXP R_x_min, SEXP R_y_min, SEXP R_x_max, SEXP R_y_max){
	// grab vector length
	int n = INTEGER(coerceVector(R_n, INTSXP))[0];
  	double bw = REAL(coerceVector(R_bw, REALSXP))[0];
  	bool torus = LOGICAL(coerceVector(R_torus, LGLSXP))[0];
  	double x_min = REAL(coerceVector(R_x_min, REALSXP))[0];
  	double y_min = REAL(coerceVector(R_y_min, REALSXP))[0];
	double x_max = REAL(coerceVector(R_x_max, REALSXP))[0];
  	double y_max = REAL(coerceVector(R_y_max, REALSXP))[0];
	// digest R_X and R_Y
	R_X=coerceVector(R_X, REALSXP);
	R_Y=coerceVector(R_Y, REALSXP);
	// output vector of the sums of toroidal distances for each individual
	SEXP outmat; PROTECT(outmat=allocVector(REALSXP, n));
	// assign pointers for the R_X, R_Y, and outmat vectors
	double *X, *Y, *out;
	X = REAL(R_X);
	Y = REAL(R_Y);
	out = REAL(outmat);
	// initialise variables to store the toroidal distance, distTemp;
	// and the density around an individual, w.
	double distTemp, w;
	// for each individual ii calculate the sum of toroidal distances between all individuals jj from 0 to n
	for (int ii=0; ii<n; ii++){
		out[ii] = 0; // initialise density at zero to prepare for adding densities for all pairwise densities
		for (int jj=0; jj<n; jj++){
			if (torus){
				distTemp = ToroidalDistance(X[ii], Y[ii], X[jj], Y[jj], x_min, y_min, x_max, y_max); // calculate toroidal distance
			} else {
				distTemp = NonToroidalDistance(X[ii], Y[ii], X[jj], Y[jj]); // calculate non-toroidal distance
			}
			if (y_min == y_max){
				w = norm1D(distTemp, bw); // calculate the density as a transformation of toroidal distance into the distance from the mean of a univariate normal distribution and calculating the corresponding probability density.
			} else {
				w = norm2D(distTemp, bw); // calculate the density as a transformation of toroidal distance into the distance from the mean of a bivariate normal distribution and calculating the corresponding probability density.
			}
			out[ii] += w; // the desnsity around individual ii is the sum of all the densities given the pairwise toroidal distances between individual ii and all other individuals
		}
	}
	UNPROTECT(1); // what does this do?
	return(outmat);
}

// Given a popmatrix, finds a mate for each individual, 2D case
// Identify the male mate of each female, 
//		if nM > nF then output is of length nF with count of -99 <= nF
//		if nM < nF then output is of length nF with count of -99 <= nM
//		also if nM < nF, allows for the same male to mate with multiple females
// NOTE: for 1D simulations set y_min == y_max
SEXP mate (SEXP R_fX, SEXP R_fY, SEXP R_mX, SEXP R_mY, SEXP R_nF, SEXP R_nM, SEXP R_bw, SEXP R_torus, SEXP R_x_min, SEXP R_y_min, SEXP R_x_max, SEXP R_y_max){
	R_fX=coerceVector(R_fX, REALSXP);
	R_fY=coerceVector(R_fY, REALSXP);
	R_mX=coerceVector(R_mX, REALSXP);
	R_mY=coerceVector(R_mY, REALSXP);
	double bw = REAL(coerceVector(R_bw, REALSXP))[0];
	int nF = INTEGER(coerceVector(R_nF, INTSXP))[0];
	int nM = INTEGER(coerceVector(R_nM, INTSXP))[0];
  	bool torus = LOGICAL(coerceVector(R_torus, LGLSXP))[0];
	double x_min = REAL(coerceVector(R_x_min, REALSXP))[0];
  	double y_min = REAL(coerceVector(R_y_min, REALSXP))[0];
	double x_max = REAL(coerceVector(R_x_max, REALSXP))[0];
  	double y_max = REAL(coerceVector(R_y_max, REALSXP))[0];

	SEXP mate; PROTECT(mate=allocVector(INTSXP, nF));
	
	const double small_num = 1e-8;
	double *fX, *fY, *mX, *mY;
	int *M;
	fX = REAL(R_fX);
	fY = REAL(R_fY);
	mX = REAL(R_mX);
	mY = REAL(R_mY);
	M = INTEGER(mate);
	double w, rnum, Dtemp, d, bw_factor;
	//double distTemp[nF][nM];// = 
	
	GetRNGstate();
	
	// find a mate for each female by randomly choosing from all the males 3*bw for for 2D and 1D (used to be within 14*bw in 1D and 3*bw in 2D.)
	bw_factor = 3;
	// if (y_min == y_max){
	// 	bw_factor = 3; /// SET TO 3 testing 2021/11/19
	// } else {
	// 	bw_factor = 3;
	// }
	for (int ii=0; ii<nF; ii++){ // step through females
		Dtemp = 0; //D[ii]-norm(0, bw);	
		for (int jj=0; jj<nM; jj++){ //step through males
			if (torus){
				// calculate distances on flat torus
				d = ToroidalDistance(fX[ii], fY[ii], mX[jj], mY[jj], x_min, y_min, x_max, y_max);
			} else {
				d = NonToroidalDistance(fX[ii], fY[ii], mX[jj], mY[jj]);
			}
			// exclude mates at threshold distance 14*bw in 1D and 3*bw in 2D
			if (d > bw_factor*bw) continue;
			// add nearby males to total male density
			if (y_min == y_max){
				Dtemp += norm1D(d, bw);
			} else {
				Dtemp += norm2D(d, bw);
			}
		}
		if (Dtemp==0) {M[ii]=-99; continue;} // catch the lonely ones
	
		rnum = unif_rand();// Random number between 0-1
		w=0;
		int jj = 0;
		do
		{  //finds the R vector index of the individual's mate
			if (torus){
				// calculate distances on flat torus
				d = ToroidalDistance(fX[ii], fY[ii], mX[jj], mY[jj], x_min, y_min, x_max, y_max);
			} else {
				d = NonToroidalDistance(fX[ii], fY[ii], mX[jj], mY[jj]);
			}
		  	// exclude mates at threshold distance 14*bw in 1D and 3*bw in 2D
			if (d > bw_factor*bw) {jj++; continue;}
		  	//cumulant, should add to one eventually
			if (y_min == y_max){
				w += norm1D(d, bw)/Dtemp;
			} else {
				w += norm2D(d, bw)/Dtemp;
			}
			//printf("ind=%i rnum=%e cumulant=%e lessthan=%i catch=%i\n", ii, rnum, w, w<rnum, (rnum-w)>small_num);
			jj++;
		} while (w<rnum && (rnum-w)>small_num); // Can we get away with just (rnum-w) > small_num? Since initially w is always <= to rnum, i.e w=0 and rnum=(0, 1)

		M[ii] = jj;	
	}
	PutRNGstate();
	UNPROTECT(1);
	return(mate);
}
