#include "main.hh"

/* This code represents a C++ encoded version of the STASH model first
 * published by Martin Sykes et al. (1992). The model operates as per the
 * formulations specified in that paper, but has been rewritten in C++ so that
 * it can be manipulated through the R Statistical language via the RCPP
 * library. This retains the speed of the original code and allows model
 * outputs to be directly analysed.
 *
 * Grided climate data is passed from R as SEXP objects and then accessed
 * 'like' a vector<type> through the RCPP syntactic sugar, such that model runs
 * as if from a direct binary. This requires that all passed SEXP objects from
 * R must be correctly formatted data frame or matrix types.
 *
 * For each temperature, rainfall and sunlight hour matrices passed and
 * translated into NumericMatrix types, the following structure holds:
 * => cols 01:02 are latitude and longitude
 * => cols 03:15 are monthly averages values of the meteorology
 *
 * GridCell objects contain private members storing the various time-dependent
 * measure of the three drivers plus the derived climatologies that have been
 * evaluated through an energy bucket model accessed from the Energy class.
 */

RcppExport SEXP gridStash( const SEXP R_gtc, const SEXP R_gpr, const SEXP R_gfs, const SEXP R_gcChar ) {

    NumericMatrix   temp(R_gtc), prec(R_gpr), fsol(R_gfs), gcChar(R_gcChar);

    NumericVector   lon = fsol(_,0), lat = fsol(_,1), elev = gcChar(_,2), fcap = gcChar(_,3);

    unsigned int    mn = 12,
                    ncol = mn + 2;
    unsigned long   ll, ncell = fsol.nrow();

    NumericMatrix   gTOT    (ncell,nvar+2),
                    gAET    (ncell,14), gEET    (ncell,14), gPET    (ncell,14),
                    gDET    (ncell,14), gPAR    (ncell,14), gMI     (ncell,14),
                    gALPHA  (ncell,14), gGDD0   (ncell,14), gGDD5   (ncell,14),
                    gGDD10  (ncell,14), gCHILL  (ncell,14), gRO     (ncell,14);


    vector< vector<float> >   tair, rain, fsun;
    GridCell        gridCell;

    // 2D storage vectors for holding values stored in R objects
    tair.resize( ncol, vector<float>( ncell, 0 ) );
    rain.resize( ncol, vector<float>( ncell, 0 ) );
    fsun.resize( ncol, vector<float>( ncell, 0 ) );
    // run through each grid cell and perform daily time-series calculations
    for( ll=0; ll<ncell; ll++ ) {

        gridCell.initialise_State( ll, lat(ll), lon(ll), elev(ll), fcap(ll) );
        gridCell.initialise_Force( tair[ll], rain[ll], fsun[ll] );

        stash_Model( gridCell )

        // convert back to SEXP objects for export to R
        assign_R_total( gridCell, ll, gTOT );
        assign_R_month( gridCell, ll, gAET,     &GridCell::get_month_AET     );
        assign_R_month( gridCell, ll, gEET,     &GridCell::get_month_EET     );
        assign_R_month( gridCell, ll, gPET,     &GridCell::get_month_PET     );
        assign_R_month( gridCell, ll, gDET,     &GridCell::get_month_DET     );
        assign_R_month( gridCell, ll, gPAR,     &GridCell::get_month_PAR     );
        assign_R_month( gridCell, ll, gRO,      &GridCell::get_month_RUN     );
        assign_R_month( gridCell, ll, gMI,      &GridCell::get_month_MI      );
        assign_R_month( gridCell, ll, gALPHA,   &GridCell::get_month_ALPHA   );
        assign_R_month( gridCell, ll, gGDD0,    &GridCell::get_month_GDD0    );
        assign_R_month( gridCell, ll, gGDD5,    &GridCell::get_month_GDD5    );
        assign_R_month( gridCell, ll, gGDD10,   &GridCell::get_month_GDD10   );
        assign_R_month( gridCell, ll, gCHILL,   &GridCell::get_month_CHILL   );
    }
    // all matrices are passed back to R as a list type
    return List::create(
            _("Total")  = gTOT, _("AET")    = gAET, _("EET")    = gEET,     _("PET")    = gPET,      _("DET") = gDET,
            _("PAR")    = gPAR, _("MI")     = gMI,  _("RO")     = gRO,      _("Alpha")  = gALPHA,
            _("GDD0")   = gGDD0,_("GDD5")   = gGDD5,_("GDD10")  = gGDD10,   _("Chill")  = gCHILL
            );
}

void stash_Model( GridCell &gc, float tair, float rain, float fsun ) {
    // create program objects
    Energy  energy;
    bool    check_pass;
    // check to see that the grid cell has no missing values (must have a complete year's worth of data)
    check_pass = missing_Value_Check( tair, rain, fsun );
    // if TRUE then run normally
    if( check_pass ) {
        // store monthly climate drivers in the gridCell object
        insert_Meteorologies( gridCell, tair, rain, fsun )
        // derive daily measurements from monthly values using linear interpolation
        interpolate_Daily( gridCell );
        // store pass value
        gridCell.set_isNull_Flag( check_pass );
        // do water balance calculations
        energy.waterBucket( gridCell );
        // perform monthly and annual sums
        gridCell.calculate_Climatologies();
        // echo to user
        cout << "STASHING >> Grid Cell :" << gridCell.get_Cell() << endl;
    // if FALSE then set all calculations for this cell to -999 (ALMA format standard)
    } else {
        gridCell.set_isNull_Flag( check_pass );
        // assign missing values to all climatologies
        gridCell.set_Missing_Value();
        // echo to user
        cout << "NOT STASHING >> Null Cell :" << gridCell.get_Cell() << endl;
    }
}

void interpolate_Daily( GridCell &gc ) {
// linearly interpolate monthly values to daily values for all climate drivers
    gc.linear_interp( gc, gc.get_month_FSUN(),  &GridCell::set_day_FSUN );
    gc.linear_interp( gc, gc.get_month_TEMP(),  &GridCell::set_day_TEMP );
    gc.linear_interp( gc, gc.get_month_PPT (),  &GridCell::set_day_PPT  );
}



/*
 * Overload functions to make semi-polymorphic assignment to SEXP objects that can be sent back to R
 */

void assign_R_total( GridCell &gc, const unsigned long ll, NumericMatrix yGrid ) {
// store all yearly climatologies into a numeric vector that can be passed back to R
    NumericVector   year_vars(nvar+2);
    year_vars = NumericVector::create(
            gc.get_Lon(),   gc.get_Lat(),
            gc.get_year_AET(),  gc.get_year_EET(),  gc.get_year_PET(),
            gc.get_year_DET(),  gc.get_year_PAR(),  gc.get_year_MI(),
            gc.get_year_ALPHA(),gc.get_year_TEMP(), gc.get_year_PPT(),
            gc.get_year_FSUN(), gc.get_year_RUN(),
            gc.get_year_GDD0(),
            gc.get_year_GDD5(),
            gc.get_year_GDD10(),
            gc.get_year_CHILL()
            );
    for( unsigned int i=0; i<nvar+2; i++ ) {
        yGrid(ll,i) = year_vars(i);
    }
}

void assign_R_month( GridCell &gc, const unsigned long ll, NumericMatrix mGrid, float (GridCell::*fget)(const int) const ) {
// create a new matrix for the climatology by passing grid cell coordinate and real value
    mGrid(ll,0) = gc.get_Lon();
    mGrid(ll,1) = gc.get_Lat();
    for( int i=0; i<gc.get_MLEN(); i++ ) {
        mGrid(ll,i+2) = (gc.*fget)(i);
    }
}

void assign_R_month( GridCell &gc, const unsigned long ll, NumericMatrix mGrid, int (GridCell::*fget)(const int) const ) {
// create a new matrix for the climatology by passing grid cell coordinate and integer value
    mGrid(ll,0) = gc.get_Lon();
    mGrid(ll,1) = gc.get_Lat();
    for( int i=0; i<gc.get_MLEN(); i++ ) {
        mGrid(ll,i+2) = (gc.*fget)(i);
    }
}

