#include "energy.hh"

RcppExport SEXP sunshine( const SEXP r_sun ) {
// Calculates the total radiation at the top of the atmosphere based on
// latitude and day of the year and the consequent monthly sunshine hours.
// -- used for gridded datasets
    NumericMatrix   swr(r_sun);

    const int   mdays[] = {31,28,31,30,31,30,31,31,30,31,30,31},
                nmon    = 12;

    int     dn, ll, i, j, ncell = swr.nrow();
    float   rtoa;

    // create energy object
    Energy          energy;

    NumericMatrix   mrtoa( ncell, nmon+2 );
    NumericMatrix   mrsun( ncell, nmon+2 );
    // fill matrices with zeroes
    fill( mrtoa.begin(), mrtoa.end(), 0 );
    fill( mrsun.begin(), mrsun.end(), 0 );
    // allocate lat and lon columns
    mrtoa(_,0) = swr(_,0);
    mrtoa(_,1) = swr(_,1);
    mrsun(_,0) = swr(_,0);
    mrsun(_,1) = swr(_,1);

    // do radiation calculations
    for( ll=0, dn=0; ll<ncell; ll++ ) {
        cout << "SHINING >> Grid Cell :" << ll+1 << endl;
        for( i=0; i<nmon; i++ ) {
            if( swr(ll,i+2)!=-999. ) {
                for( j=0; j<mdays[i]; j++, ++dn ) {
                    rtoa = energy.terrestrial_radiation( dn, swr(ll,1) );     // column index 1 of the swr matrix is the lat column
                    mrtoa(ll,i+2) += rtoa/mdays[i];
                }
                mrsun(ll,i+2) = ( (swr(ll,i+2)/mrtoa(ll,i+2))-0.25 )/0.5;
            } else {
                mrtoa(ll,i+2) = -999.;
                mrsun(ll,i+2) = -999.;
            } // end day
        } // end month
    } // end cell
    return List::create( _("SunHr")=mrsun, _("RadTop")=mrtoa );
} // end program

//================================================================================

RcppExport SEXP sunray( const SEXP R_swr, const SEXP R_mon, const SEXP R_lat ) {
// Calculates the total radiation at the top of the atmosphere based on
// latitude and day of the year and the consequent monthly sunshine hours.
// -- used for non-gridded datasets
    const int   mlen[] = {31,28,31,30,31,30,31,31,30,31,30,31},
                mbeg[] = {0,31,59,90,120,151,181,212,243,273,304,334};

    NumericVector   swr(R_swr), mon(R_mon), lat(R_lat);
    int             ix, dn,
                    vlen = swr.size();
    NumericVector   mrsun(vlen), mrtoa(vlen);
    // fill matrices with zeroes
    fill( mrtoa.begin(), mrtoa.end(), 0 );
    fill( mrsun.begin(), mrsun.end(), 0 );

    // create energy object
    Energy          energy;

    for( int i=0; i<vlen; i++ ) {
        ix  = mon(i)-1;
        dn  = mbeg[ix]+1;
        if( swr(i)!=-999. ) {
            for( int j=0; j<mlen[ix]; j++, ++dn ) {
                mrtoa(i) += energy.terrestrial_radiation( dn, lat(i) )/mlen[ix];
            }
            mrsun(i) = ( swr(i)/mrtoa(i)-0.25 )/0.5;
        } else {
            mrtoa(i)    = -999.;
            mrsun(i)    = -999.;
        } // end day
    }
    return List::create( _("SunHr")=mrsun, _("RadTop")=mrtoa );
} // end sunray

//================================================================================

RcppExport SEXP pointStash( const SEXP R_tc, const SEXP R_pr, const SEXP R_fs, const SEXP R_mon, const SEXP R_lat, const SEXP R_fcap, const SEXP R_elev ) {
// takes climate data at the point measure and determines bioclimatic information
// no monthly interpolation nor water balance is performed - this is purely determined
// monthly totals and averages
// ** assumes that there is no limitation to soil water supply
    Energy          energy;

    const int       mlen[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    IntegerVector   mon(R_mon);
    NumericVector   temp(R_tc), prec(R_pr), fsun(R_fs),
                    lat(R_lat), fcap(R_fcap), elev(R_elev);
    int             day, monl,
                    N = temp.size();
    float           spl, par, aet, eet, pet, det;
    vector<float>   mpar(N), maet(N), meet(N), mpet(N), mdet(N), mmi(N), malpha(N);
    // calculate
    for( int i=0; i<N; i++ ) {
        monl        = mlen[ mon(i)-1 ];
        // use the middle of the month
        day         = monl/2;
        // soil plant limiting multiplier - if exceeds field capacity then AET is at a maximum
        spl         = min( prec(i)/fcap(i), 1.0 );
        // average energy exchange for the month is calculated
        //if( temp(i)!=-9999. || fsun(i)!=-9999. || lat(i)!=-9999. ) {
        if( temp(i)!=-999. and fsun(i)!=-999. and lat(i)!=-999. ) {
            energy.evaporate( temp(i), fsun(i), lat(i), spl, day, par, aet, eet, pet, det );
            mpar[i]     = par*monl;
            maet[i]     = aet*monl;
            meet[i]     = eet*monl;
            mpet[i]     = pet*monl;
            mdet[i]     = det*monl;
            // monthly totals determined from month length
            malpha[i]   = aet/eet;
            mmi[i]      = prec(i)/meet[i];
        } else {
            mpar[i]     = -999.;
            maet[i]     = -999.;
            meet[i]     = -999.;
            mpet[i]     = -999.;
            mdet[i]     = -999.;
            mmi[i]      = -999.;
            malpha[i]   = -999.;
        }
    }
    return List::create(
            _("AET")    = maet, _("EET")    = meet, _("PET")    = mpet,      _("DET") = mdet,
            _("PAR")    = mpar, _("MI")     = mmi,  _("Alpha")  = malpha
            );
}


