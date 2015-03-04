#include "gridcell.hh"

const float GridCell::miss_val = -999.0;

// Setters
void GridCell::set_Elevation        ( const float val ) { elev = val; }
void GridCell::set_Field_Capacity   ( const float val ) { fcap = val; }
void GridCell::set_Coordinates      ( const float rlat, const float rlon )
                                                        { lat = rlat; lon = rlon; }
void GridCell::set_Cell_Number      ( const unsigned int n )
                                                        { cell = n; }
void GridCell::set_isNull_Flag      ( const bool flag ) { isnull = flag; }
void GridCell::set_WaterBalance_Flag()                  { wbss = false; }

// Getters
float GridCell::get_Cell_Number     () const { return cell; }
float GridCell::get_Latitude        () const { return lat; }
float GridCell::get_Longitude       () const { return lon; }
float GridCell::get_Field_Capacity  () const { return fcap; }
float GridCell::get_Elevation       () const { return elev; }

// FUNCTIONS
void GridCell::initialise_Cell( const unsigned int cell_num, const float lat, const float lon, const float elev, const float fcap, const bool isleap ) {
// set the cell number and other geo-referential site characteristics
    set_Cell_Number    ( cell_num );
    set_Coordinates    ( lat, lon );
    set_Elevation      ( elev );
    set_Field_Capacity ( fcap );
    set_isLeap_Flag    ( isleap );
}

void GridCell::initialise_State() {
// initialise the state at time t=0 for the model. Vector lengths are
// initialised and all values set to 0
    if( get_LEAP()==true ) {
        init_Leap();
    } else {
        init_noLeap();
    }
    init_Month();
    resIn_Year(0.0);
    // first day water content is set to the field capacity
    set_day_SMC(fcap,0);
}

bool GridCell::initialise_Force( const vector<float> tair, const vector<float> rain, const vector<float> fsun ) {
// initialise the meteorologies for the gridCell that will drive the model.
// Data must pass a missing value check in order for the model to be run on the
// cell. If values are missing then set forcing to -999.0 and return false to
// skip the cell from running in the model
    bool pass;

    for( unsigned int mn=0; mn<get_MLEN(); mn++ ) {
        set_month_TEMP( tair[mn], mn );
        set_month_PPT ( rain[mn], mn );
        set_month_FSUN( fsun[mn], mn );
    }
    // DO MISSING VALUE CHECK HERE!
    pass = missing_ValueCheck_Flag();
    if( pass==true ) {
        return true;
    } else {
        set_Missing_Value();
        return false;
    }
}

void GridCell::set_Missing_Value() {
// the above missing value check has failed, therefore initialise all values to
// the missing value
    reset_Day   ( miss_val );
    reset_Month ( miss_val );
    resIn_Year  ( miss_val );
}


bool GridCell::missing_ValueCheck_Flag() {
// search through each of the meteorologies looking for missing values (-999),
    vector<float>::iterator check_tair, check_rain, check_fsun;
    // look for missing values for all three drivers individually
    check_tair = find( get_month_TEMP().begin(), get_month_TEMP().end(), miss_val );
    check_rain = find( get_month_FSUN().begin(), get_month_FSUN().end(), miss_val );
    check_fsun = find( get_month_PPT ().begin(), get_month_PPT ().end(), miss_val );
    // if there are no missing values then check passes
    if( (*check_tair != miss_val)  | (*check_rain != miss_val) | (*check_fsun != miss_val) ) {
        return true;
    // else the check fails
    } else {
        return false;
    }
}

// linear interpolation -- convert monthly to daily
void GridCell::linear_interp( GridCell &gc, const vector<float> mly12, void (GridCell::*fset)(const float, const int) ) {
    // Vector of days in month (includes bookend months - Dec year before and Jan year after)
    const float m14_noleap[] = {31,31,28,31,30,31,30,31,31,30,31,30,31,31},
                m14_leap[] = {31,31,29,31,30,31,30,31,31,30,31,30,31,31},
                *m14days;

    // Matrix of 12 month columns and N grid cell rows
    vector<float>   mly14;
    float           mbef, maft, inc, d15;
    int             i, j, dn, dbef, daft,
                    nm12 = mly12.size(),
                    nm14 = nm12+2;

    // handle leap years
    if( get_LEAP() ) {
        m14days = m14_leap;
    } else {
        m14days = m14_noleap;
    }
    // must include months either side of a year
    mly14.resize(nm14);
    mly14[0]  = mly12[11];
    mly14[13] = mly12[1];
    for( i=0; i<nm12; i++ ) {
        mly14[i+1] = mly12[i];
    }
    // Perform a linear interpolation between monthly 'mid-point' values
    for( dn=0, i=1; i<(nm14-1); i++ ) {
        for( j=0; j<m14days[i]; j++, ++dn ) {
            // first half of month
            if( j<(m14days[i]/2) ) {
                mbef = mly14[i-1];
                maft = mly14[i];
                dbef = m14days[i-1];
                daft = m14days[i];
            // second half of month
            } else {
                mbef = mly14[i];
                maft = mly14[i+1];
                dbef = m14days[i];
                daft = m14days[i+1];
            }
            inc = (maft-mbef)/((dbef+daft)/2);
            if( inc!=0 ) {
                if( j<(m14days[i]/2) ) {
                    d15 = j+(dbef/2);
                } else {
                    d15 = j-(daft/2);
                }
                (gc.*fset)( mbef+(d15*inc), dn );
            } else {
                (gc.*fset)( mbef, dn );
            }
        }
    }
}
// end routine

// ifelse hack
float GridCell::calcGDD( const float gtemp, const int gdays ) {
    if( gdays>0 ) {
        return gtemp/(float)(gdays);
    } else {
        return 0.0;
    }
}

// growing degree days
void GridCell::growDegDay() {
    for( unsigned int i=0; i<get_MLEN(); i++ ) {
//        mgdd0[i]    = calcGDD( gdd0[i], gday0[i] );
//        mgdd5[i]    = calcGDD( gdd5[i], gday5[i] );
//        mgdd10[i]   = calcGDD( gdd10[i], gday10[i] );
    }
}

// monthly sums of energy/water-use
void GridCell::monthlySums() {
//    int         dn = 0;
    for( unsigned int i=0; i<get_MLEN(); i++ ) {
        for( unsigned int j=0; j<get_MDAY(i); j++ ) {
//            mpet[i] += dpet[dn];
//            maet[i] += daet[dn];
//            meet[i] += deet[dn];
//            mdet[i] += ddet[dn];
//            mpar[i] += dpar[dn];
//            mro[i]  += dro[dn];
//            dn++;
        }
    }
}

void GridCell::annualSums() {
    // sums
//    yaet    = accumulate( maet  .begin(),   maet    .end(), 0.0 );
//    yeet    = accumulate( meet  .begin(),   meet    .end(), 0.0 );
//    ypet    = accumulate( mpet  .begin(),   mpet    .end(), 0.0 );
//    ydet    = accumulate( mdet  .begin(),   mdet    .end(), 0.0 );
//    ypar    = accumulate( mpar  .begin(),   mpar    .end(), 0.0 );
//    yro     = accumulate( mro   .begin(),   mro     .end(), 0.0 );
//    ypr     = accumulate( mpr   .begin(),   mpr     .end(), 0.0 );
//    ychill  = accumulate( mchill.begin(),   mchill  .end(), 0.0 );
//    ygdd0   = accumulate( mgdd0 .begin(),   mgdd0   .end(), 0.0 );
//    ygdd5   = accumulate( mgdd5 .begin(),   mgdd5   .end(), 0.0 );
//    ygdd10  = accumulate( mgdd10.begin(),   mgdd10  .end(), 0.0 );
//    // means
//    ytc     = accumulate( mtc   .begin(),   mtc     .end(), 0.0 ) / mlen;
//    yfs     = accumulate( mfs   .begin(),   mfs     .end(), 0.0 ) / mlen;
//    // indices
//    ymi     = ypr/yeet;
//    yalpha  = yaet/yeet;
}

// monthly indices of wet and dry
void GridCell::monthlyIndex() {
    for( unsigned int i=0; i<get_MLEN(); i++ ) {
//        mmi[i]    = mpr[i]/meet[i];
//        malpha[i] = maet[i]/meet[i];
    }
}

void GridCell::calculate_Climatologies() {
// bulk calculation of all climatologies
//    growDegDay();
//    monthlySums();
//    monthlyIndex();
//    annualSums();
}

///////////////////////////////////////////////////////////////////////////////////
///// XU & HUTCHINSON INDICES FUNCTIONS ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

////Calculate tc with tmin and tmax
//void GridCell::calc_tc() {
//    for (int i=0; i<mlen; i++){
//        mtc[i] = (mtn[i] + mtx[i])/2 ;
//    }
//}
//
//void GridCell::calc_dpr( GridCell &gc ) {
//    int dn=0, i, j;
//    float ppt;
//    const int   mdays[] = {31,28,31,30,31,30,31,31,30,31,30,31};
//    for( dn=0, i=0; i<gc.get_MLEN(); i++ ) {
//        for( j=0; j<mdays[i]; j++, dn++ ) {
//            // rainfall needs to be a fraction of the monthly total (should fix this initially)
//            ppt = gc.get_day_PPT(dn)/mdays[i];
//            gc.set_day_PPT(ppt, dn);
//        }
//    }
//}
//
//// Annual mean/sum (1 & 12 & 20)
//void GridCell::annual_Variables () {
//    ybio01  = accumulate( mtc   .begin(),   mtc     .end(), 0.0 ) / mlen;
//    ybio12  = accumulate( mpr   .begin(),   mpr     .end(), 0.0 );
//    ybio20  = accumulate( mrd   .begin(),   mrd     .end(), 0.0 ) / mlen;
//}
//
//float GridCell::max_element (vector <float> A){
//    float n = A[0];
//    for (unsigned int i=1; i<A.size(); i++){
//        if (A[i]>n) {n=A[i];}
//    }
//    return n;
//}
//
//float GridCell::min_element (vector <float> A){
//    float n = A[0];
//    for (unsigned int i=1; i<A.size(); i++){
//        if (A[i]<n) {n=A[i];}
//    }
//    return n;
//}
//
//// Quarterly variables (8,9,10,11 & 16,17,18,19 & 24,25,26,27)
//void GridCell::quarter_Variables () {
//
//    // vector of  monthly pr/tc/fs (include Nov Dec before and Jan Feb after)
//    vector <float> mpr16(16), mtc16(16), mrd16(16);
//    mpr16[0] = mpr[10]; mpr16[1] = mpr[11]; mpr16[14] = mpr[0]; mpr16[15] = mpr[1];
//    mtc16[0] = mtc[10]; mtc16[1] = mtc[11]; mtc16[14] = mtc[0]; mtc16[15] = mtc[1];
//    mrd16[0] = mrd[10]; mrd16[1] = mrd[11]; mrd16[14] = mrd[0]; mrd16[15] = mrd[1];
//    for( unsigned int i=0; i<mpr.size(); i++ ) {
//        mpr16[i+2] = mpr[i];
//        mtc16[i+2] = mtc[i];
//        mrd16[i+2] = mrd[i];
//    }
//
//    vector <float> qtrPr(14); // sum of precipitation for all possible quarters of the year
//    vector <float> qtrTc(14); // mean temperature for all possible quarters of the year
//    vector <float> qtrRd(14); // mean radiation for all possible quarters of the year
//    for (unsigned int i=0; i<qtrPr.size(); i++) {
//        qtrPr[i]=mpr16[i]+mpr16[i+1]+mpr16[i+2];
//        qtrTc[i]=(mtc16[i]+mtc16[i+1]+mtc16[i+2])/3;
//        qtrRd[i]=(mrd16[i]+mrd16[i+1]+mrd16[i+2])/3;
//    }
//
//    // max and min value of qtrPr and qtrT
//    float maxPr = max_element(qtrPr);
//    float minPr = min_element(qtrPr);
//    float maxTc = max_element(qtrTc);
//    float minTc = min_element(qtrTc);
//
//    // first month of the wettest/driest/warmest/coldest quarter
//    int wet=0; while (maxPr!=qtrPr[wet]){wet++;}
//    int dry=0; while (minPr!=qtrPr[dry]){dry++;}
//    int warm=0; while (maxTc!=qtrTc[warm]){warm++;}
//    int cold=0; while (minTc!=qtrTc[cold]){cold++;}
//
//    // Temperature variables
//    ybio08 = qtrTc[wet]; // BIO08 Mean Temperature of Wettest Quarter
//    ybio09 = qtrTc[dry]; // BIO09 Mean Temperature of Driest Quarter
//    ybio10 = qtrTc[warm]; // BIO10 Mean Temperature of Warmest Quarter
//    ybio11 = qtrTc[cold]; // BIO11 Mean Temperature of Coldest Quarter
//
//    // Precipitation variables
//    ybio16 = qtrPr[wet]; // BIO16 Precipitation of Wettest Quarter
//    ybio17 = qtrPr[dry]; // BIO17 Precipitation of Driest Quarter
//    ybio18 = qtrPr[warm]; // BIO18 Precipitation of Warmest Quarter
//    ybio19 = qtrPr[cold]; // BIO19 Precipitation of Coldest Quarter
//
//    // Radiation variables
//    ybio24 = qtrRd[wet] ;// BIO24 Radiation of Wettest Quarter
//    ybio25 = qtrRd[dry] ;// BIO25 Radiation of Driest Quarter
//    ybio26 = qtrRd[warm] ;// BIO26 Radiation of Warmest Quarter
//    ybio27 = qtrRd[cold] ;// BIO27 Radiation of Coldest Quarter
//}
//
//// Period variables (5,6,7 & 13,14 & 21,22)
//void GridCell::period_Variables ( GridCell &gc) {
//    // max and min value of mpr and mtc
//    float max_mpr = max_element(mpr);
//    float min_mpr = min_element(mpr);
//    float max_mtc = max_element(mtc);
//    float min_mtc = min_element(mtc);
//    float max_mrd = max_element(mrd);
//    float min_mrd = min_element(mrd);
//
//    // month of the wettest/driest/warmest/coldest/high and low rad period
//    int wet=0;  while (max_mpr!=mpr[wet]) {wet++;}
//    int dry=0;  while (min_mpr!=mpr[dry]) {dry++;}
//    int warm=0; while (max_mtc!=mtc[warm]){warm++;}
//    int cold=0; while (min_mtc!=mtc[cold]){cold++;}
//    int high=0; while (max_mrd!=mrd[high]){high++;}
//    int low=0;  while (min_mrd!=mrd[low]) {low++;}
//
//    // Temperature variables
//    ybio05 = gc.get_month_TMAX(warm) ;// BIO05 Maximum Temperature of Warmest Period
//    ybio06 = gc.get_month_TMIN(cold) ;// BIO06 Minimum Temperature of Coldest Period
//
//    ybio07 = ybio05 - ybio06 ; // BIO07 Temperature Annual Range (BIO05 - BIO06)
//
//    // Precipitation variables
//    ybio13 = mpr[wet] ;// BIO13 Precipitation of Wettest Period
//    ybio14 = mpr[dry] ;// BIO14 Precipitation of Driest Period
//
//    // Radiation variables
//    ybio21 = mrd[high] ;// BIO21 Highest Period Radiation
//    ybio22 = mrd[low] ;// BIO22 Lowest Period Radiation
//}
//
//// Seasonality variables (standard deviation * 100) (4,15,23)
//void GridCell::seasonal_Variables () {
//    const int   mdays[] = {31,28,31,30,31,30,31,31,30,31,30,31};
//    int dn=0;
//    for (int i=0; i<mlen; i++){
//        float temp=0, prec=0, rad=0, pmean1=(mpr[i]+mdays[i])/mdays[i];
//        for (int j=0; j<mdays[i]; j++){
//            temp += (mtc[i] - dtc[dn])*(mtc[i] - dtc[dn]);
//            prec += (pmean1 - dpr[dn])*(pmean1 - dpr[dn]);
//            rad  += (mrd[i] - drd[dn])*(mrd[i] - drd[dn]);
//            dn++;
//        }
//        mbio04[i] = sqrt(temp/(mdays[i]-1))*100 ;// BIO04 Temperature Seasonality
//        mbio15[i] = (sqrt(prec/(mdays[i]-1))/pmean1)*100 ;// BIO15 Precipitation Seasonality
//        mbio23[i] = sqrt(rad/(mdays[i]-1))*100 ;// BIO23 Radiation Seasonality
//    }
//
//    float temp2=0, prec2=0, rad2=0, pmean=(ybio12+mlen)/mlen;
//    for (int k=0; k<mlen; k++){
//        temp2 += (ybio01 - mtc[k])*(ybio01 - mtc[k]);
//        prec2 += (pmean - (mpr[k]+1))*(pmean - (mpr[k]+1));
//        rad2  += (ybio20 - mrd[k])*(ybio20 - mrd[k]);
//    }
//    ybio04 = sqrt(temp2/(mlen-1))*100 ;// BIO04 Temperature Seasonality
//    ybio15 = (sqrt(prec2/(mlen-1))/pmean)*100 ;// BIO15 Precipitation Seasonality
//    ybio23 = sqrt(rad2/(mlen-1))*100 ;// BIO23 Radiation Seasonality
//}
//
//// Mean diurnal Range (2)
//void GridCell::mean_Diurnal_Range ( GridCell &gc){
//    float tmin, tmax, b2=0;
//    for (int i=0; i<mlen; i++) {
//        // getting
//        tmin    = gc.get_month_TMIN(i);
//        tmax    = gc.get_month_TMAX(i);
//        // calculations
//        mbio02[i] = tmax-tmin;
//        b2 = b2 +mbio02[i];
//    }
//    ybio02 = b2 / mlen;
//}
//
//// Isothermality ((bio02/bio07)*100) (3)
//void GridCell::isothermality (){
//    for (unsigned int i=0; i<mbio03.size();i++){
//        mbio03[i] = (mbio02[i]/ybio07)*100;
//    }
//    ybio03 = (ybio02/ybio07)*100;
//}
//
