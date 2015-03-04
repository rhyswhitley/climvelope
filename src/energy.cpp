// functions to determine monthly evaporation products inclusive of potential, actual evapotransM_PIration, aridity and
// moisture indices and soil moisture controls.

#include "energy.hh"

float Energy::solar_const( const int day ) {
// calculates the solar constant based on day
    static const float eccen = 0.01675;
    float arg, o_solc;
    arg = solar_angle(day);
    o_solc = solc*(1.+2.*eccen*cos(arg*d2r));
    return o_solc;
}

float Energy::change_solar_angle( const int day ) {
// calculates the change in solar angle (solar declination)
    float arg, o_delta;
    arg = solar_angle(day+10);
    o_delta = -23.4*cos(arg*d2r);
    return o_delta;
}

float Energy::slope_vapor_temperature( const float tair ) {
// calculates the slope of the water vaporization curve
    static const float k1 = 610.78, k2 = 17.269, k3 = 273.3;
    float denom, esat, slope;
    denom   = k3+tair;
    esat    = k1*exp(k2*tair/denom);
    slope   = k2*k3*esat/pow(denom,2);
    return slope;
}

float Energy::calc_parameter_hs( float const x, float const y, float const lat, float const delta ) {
    float Hs;
    if( x>y ) {
        Hs = M_PI;
    } else if( x< -1*y ) {
        Hs = 0.0;
    } else {
        Hs = acos(-1.*tan(lat*d2r)*tan(delta*d2r));
    }
    return Hs;
}

void Energy::calc_parameters_h01( const float u, const float v, const float spl, float& h0, float& h1 ) {
    float arg;
    if( (spl-u)>=v ) {
//  supply exceeds demand
        h1 = 0.;
    } else if( (spl-u)<= -1*v ) {
//  demand exceeds supply
        h0 = M_PI;
        h1 = M_PI;
    } else {
        arg = (spl-u)/v;
        h1 = acos(arg);
    }
    if( u>=v ) {
//  polar day
        h0 = M_PI;
    } else if( u<= -1*v ) {
//  polar night
        h0 = 0.0;
        h1 = 0.0;
    } else {
//  normal day and night
        arg = -(u/v);
        h0 = acos(arg);
    }
}

void Energy::energy_lookup_table( const float temp, float& gamma, float& lambda ) {
// reference table for temperature standardised values of gamma and lambda
    static const float tbase[] = {-5,0,5,10,15,20,25,30,35,40,45};
    static const float gbase[] = {64.6,64.9,65.2,65.6,65.9,66.1,66.5,66.8,67.2,67.5,67.8};
    static const float lbase[] = {2.531,2.501,2.489,2.477,2.465,2.454,2.422,2.43,2.418,2.406,2.394};
    bool    found   = false;
    int     i       = 0;
    while( found==false and i<12 ) {
        if( temp<=tbase[i] ) {
            gamma   = gbase[i];
            lambda  = lbase[i];
            found   = true;
        } else {
            gamma   = gbase[10];
            lambda  = lbase[10];
            found   = true;
        }
        i++;
    }
}

float Energy::terrestrial_radiation( const int day, const float lat ) {
// Calculates the total radiation at the top of the atmosphere based on latitude and day of the year.
    float   msolc, delta, x, y, hs, rtoa_MJ;
    msolc   = solar_const(day);
    delta   = change_solar_angle(day);
    // get x, y projections
    euclidean_projection( lat, delta, x, y );
    hs      = calc_parameter_hs( x, y, lat, delta );
    rtoa_MJ = par_atmos_top( msolc, hs, x, y )*nsecs*1e-6;
    return rtoa_MJ;
}

float Energy::par_atmos_top( const float msolc, const float hs, const float x, const float y) {
// calculates PAR at the top of the atmosphere (W m-2)
    float rta;
    rta = msolc*(hs*x+sin(hs)*y)/M_PI;
    return rta;
}

void Energy::euclidean_projection( const float lat, const float delta, float& x, float& y) {
// solar coordinate projections
    x = sin(lat*d2r)*sin(delta*d2r);
    y = cos(lat*d2r)*cos(delta*d2r);
}

float Energy::equilibrium_evaporaiton( const float h0, const float u, const float v) {
// equilibrium evaporation
    float pet;
    pet = (2.*u*h0+2.*v*sin(h0))/(M_PI/12.);
    return pet;
}

float Energy::actual_evaporation( const float pet, const float h1, const float spl, const float u, const float v, const float dcon ) {
// actual evaporation
    float aet;
    aet = pet + 2*(spl*h1-u*h1-v*sin(h1) )/(M_PI/12) - dcon;
    return aet;
}

float Energy::photo_active_radiation( const float rabs, const float hs, const float rsin, const float rcos ) {
// photosynthetically active radiation
    float par;
    par = solc/0.22*rabs*( hs*rsin + sin(hs)*rcos )/(2*M_PI);
    return par;
}

//void Energy::sunshine( const float swr, const float rtoa, float& fsun, float& totrad ) {
//// Calculates the total radiation at the top of the atmosphere based on
//// latitude and day of the year and the consequent monthly sunshine hours.
//// -- used for gridded datasets
//    const int   mdays[] = {31,28,31,30,31,30,31,31,30,31,30,31},
//                nmon    = 12;
//
//    int     dn, ll, i, j, ncell = swr.nrow();
//    float   rtoa;
//
//    for( j=0; j<mdays[i]; j++, ++dn ) {
//        rtoa = terrestrial_radiation( dn, lat(ll,1) );     // column index 1 of the swr matrix is the lat column
//        mrtoa(ll,i) += rtoa/mdays[i];
//    }
//
//    // do radiation calculations
//    for( ll=0, dn=0; ll<ncell; ll++ ) {
//        cout << "SHINING >> Grid Cell :" << ll+1 << endl;
//        for( i=0; i<nmon; i++ ) {
//            if( swr(ll,i)!=-999. ) {
//                fsun = sunlight_hours( swr, rtoa );
//            } else {
//                mrtoa(ll,i+2) = -999.;
//                mrsun(ll,i+2) = -999.;
//            } // end day
//        } // end month
//    } // end cell
//    return List::create( _("SunHr")=mrsun, _("RadTop")=mrtoa );
//} // end program
//

void Energy::grid_evaporate_wrapper( GridCell &gc, const float spl, const int day ) {
// calculate the exchange of energy between the earth's surface and the atmosphere
    unsigned int    dn = day-1;
    float           tair, fsun, lat, aet, pet, eet, det, par;
// getting
    tair    = gc.get_day_TEMP(dn);
    fsun    = gc.get_day_FSUN(dn); // make sure that this is already between 0 & 1
    lat     = gc.get_Lat();
// calculations
    evaporate( tair, fsun, lat, spl, day, par, aet, eet, pet, det );
// setting
    gc.set_day_AET( aet, dn );
    gc.set_day_EET( eet, dn );
    gc.set_day_PET( pet, dn );
    gc.set_day_DET( det, dn );
    gc.set_day_PAR( par, dn );
}

void Energy::evaporate( const float tair, const float fsun, const float lat, const float spl, const int day,
                float& par, float& aet, float& eet, float& pet, float& det ) {
// calculate the exchange of energy between the earth's surface and the atmosphere
    float           delta, dSat, u, v, x, y, dcon, econ,
                    msolc, rabs, rout, rsin, rcos, hs, h0, h1, gamma, lambda;
// lookup gamma and lambda values based on temperature
    energy_lookup_table( tair, gamma, lambda );
// calculations
    msolc   = solar_const(day);
    delta   = change_solar_angle(day);
    dSat    = slope_vapor_temperature(tair);
    rabs    = absorb_radiation(fsun);
    econ    = energy_convert( lambda, gamma, dSat );
// project onto Euclidean space
    euclidean_projection( lat, delta, rsin, rcos );
// energy reflected
    rout    = reflect_radiation( fsun, tair );    // long-wave radiation
// condensation
    dcon    = econ*rout;
// supply and demand for water
    x       = rabs*rsin;
    y       = rabs*rcos;
    u       = econ*(msolc*x-rout);      // net radiation
    v       = econ*msolc*y;
    hs      = calc_parameter_hs( x, y, lat, delta );
// h0 is the point in time when net radiation crosses zero
    calc_parameters_h01( u, v, spl, h0, h1 );
// calculate total daily PAR required for sum over the growing season
    par     = photo_active_radiation(rabs,hs,rsin,rcos)*1e-6*nsecs;
// potential evaporation
    eet     = equilibrium_evaporaiton(h0,u,v);
// equilibrium evaporation
    pet     = eet*1.26;
// actual evaporation
    aet     = actual_evaporation(pet,h1,spl,u,v,dcon);
// difference between equilibrium and actual
    det     = eet-aet;
}

void Energy::waterBucket( GridCell &gc ) {

    const int   mdays[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    int         i, j, dn, runc=0;
    float       ysm, lsm, fcap, tair, ppt, spl=1, smc=0, run=0;
    bool        ss = false;

    fcap    = gc.get_Fcap();
    lsm     = fcap;
    do {
        // month
        for( dn=0, i=0; i<gc.get_MLEN(); i++ ) {
            // day
            for( j=0; j<mdays[i]; j++, dn++ ) {
                if( dn==0 ) {
                    ysm = lsm;
                } else {
                    ysm = gc.get_day_SMC(dn-1);  // morning soil water equals yesterday
                }
                spl = ysm/fcap;
                grid_evaporate_wrapper( gc, spl, dn+1 );
                // rainfall needs to be a fraction of the monthly total (should fix this initially)
                ppt = gc.get_day_PPT(dn)/mdays[i];
                // soil water content is rain water minus evaporated water
                smc = ysm + ppt - gc.get_day_AET(dn);
                // determine runoff and new value
                if( smc>fcap ) {
                    run = smc - fcap;
                    gc.set_day_RUN( run, dn );
                    gc.set_day_SMC( fcap, dn );
                } else if( smc>0 & smc<fcap ) {
                    gc.set_day_RUN( 0.0, dn );
                    gc.set_day_SMC( smc, dn );
                } else {
                    gc.set_day_RUN( 0.0, dn );
                    gc.set_day_SMC( 0.0, dn );
                }
                // calculate growing degree days
                tair = gc.get_day_TEMP(dn);
                if( tair<-5.0 ) gc.set_CHILL ( i );
                if( tair>0.0  ) gc.set_GDEG0 ( tair, i );
                if( tair>5.0  ) gc.set_GDEG5 ( tair, i );
                if( tair>10.0 ) gc.set_GDEG10( tair, i );
            } // end day
        } // end month
        // checks to see if the grid's water balance is in steady-state
        if( (int)(lsm*10)!=(int)(gc.get_day_SMC(364)*10) ) {
            runc++;
            if( runc>40 ) {
                //cout << "ERROR: Water balance does not reach steady-state conditions." << endl;
                cout << "ERROR: WB:\t";
                gc.set_WaterBalance_Flag();
                break;
            }
            lsm = gc.get_day_SMC(364);
        } else {
            ss = true;
        }
    } while( ss==false );
} // end waterBucket
