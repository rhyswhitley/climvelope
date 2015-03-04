#ifndef _gridcell_H
#define _gridcell_H
#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <numeric>
#include <math.h>
#include <algorithm>
#include "energies.hh"
#include "climindex.hh"

using namespace std;

// author: Rhys Whitley
// institution: Macquarie University
// date: Jan 2015

class GridCell {
/* Objects built from this class define a geographically referenced grid cell that holds estimations
 * of climate-surface energy exchange. Energy-exchange variables are used to determine bioclimatic
 * information for a given coordinate.
 * Based off the STASH sourced code (M. Sykes 1996) and PEATSTASH source code (Gallego-Sala 2010).
 * Modified to C++ by R. Whitley.
 */
private:
    // OBJECT VARIABLES
    static const unsigned int
                    miss_val;
    static const float
                    dlen_leap,  // number of days in a leap year
                    dlen_noLeap,// number of days in a year
                    mlen,       // number of months in a year
                    ylen;       // number of years (if any)
    // Vector of days in month (includes bookend months - Dec year before and Jan year after)
    static const float
                    m14days[];

    bool            spinup, // does the model spin-up on the grid cell to equilibrate the water balance?
                    wbss,   // is the water balance of the grid cell steady-state?
                    isnull; // is the grid cell in the ocean? (therefore equals Null)
    unsigned int    cellnum;// the grid cell number
    float           lat,    // latitudinal coordinates
                    lon,    // longitudinal coordinates
                    elev,   // elevation of the grid cell
                    fcap;   // field capacity of the soil for the grid cell

public:
    // CONSTRUCTOR
    GridCell() {
        fcap    = 150.0;
        elev    = 0.0;
        spinup  = true;
        wbss    = false;
        isnull  = false;
    };

    // SETTERS
    void set_Cell_Number    ( const unsigned int n );
    void set_Coordinates    ( const float ilat,
                              const float ilon );
    void set_Elevation      ( const float value );
    void set_Field_Capacity ( const float value );
    void set_isNull_Flag    ( const bool flag );
    void set_WaterBalance_Flag
                            ( const bool flag );
    void set_SpinUp         ( const bool flag);
    void set_Missing_Value  ();

    // GETTERS
    float get_Cell_Number   () const;
    float get_Latitude      () const;
    float get_Longitude     () const;
    float get_Field_Capacity() const;
    float get_Elevation     () const;

    // FUNCTIONS
    void growDegDay         ();
    void monthlySums        ();
    void monthlyIndex       ();
    void annualSums         ();
    void calculate_Climatologies();
    void initialise_State   ();
    void initialise_Cell    ( const unsigned int cell_num, const float lat, const float lon, const float elev, const float fcap, const bool isleap );

    bool initialise_Force   ( const vector<float> tair, const vector<float> rain, const vector<float> fsun );
    bool missing_ValueCheck_Flag();

    // may move this to main
    void linear_interp      ( GridCell &gc, const vector<float> mly12, void (GridCell::*fset)(const float, const int) );

    float calcGDD( const float gtemp, const int gdays );

    // DESTRUCTOR
    ~GridCell(){};
};
#endif
