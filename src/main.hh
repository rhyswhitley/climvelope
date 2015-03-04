#ifndef _main_RCPP_MAIN_H
#define _main_RCPP_MAIN_H
#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <math.h>
#include <Rcpp.h>
#include "gridcell.hh"
#include "energy.hh"

using namespace std;
using namespace Rcpp;

static const unsigned int nvar = 15;    // the number of bioclimatic variables being returned
static const float miss_val = -999.;

RcppExport SEXP gridStash   ( const SEXP R_gtc, const SEXP R_gpr, const SEXP R_gfs, const SEXP R_gcChar );
//NumericVector   RConvert        ( vector<GridCell> &gc );
bool missing_Value_Check    ( float chk_tair, float chk_rain, float chk_fsun );
void insert_Meteorologies   ( GridCell &gc, float tair, float rain, float fsun );
void stash_Model            ( GridCell &gc, float tair, float rain, float fsun );
void interpolate_Daily      ( GridCell &gc );

// these need to eventually be in the R rwrapper file
void assign_R_total( GridCell &gc, const unsigned long ll, NumericMatrix yGrid );
void assign_R_month( GridCell &gc, const unsigned long ll, NumericMatrix mGrid, int (GridCell::*fget)(const int) const );
void assign_R_month( GridCell &gc, const unsigned long ll, NumericMatrix mGrid, float (GridCell::*fget)(const int) const );

#endif
