#ifndef _energy_H
#define _energy_H
#define _USE_MATH_DEFINES

#include <math.h>
#include "gridcell.hh"

using namespace std;
using namespace Rcpp;

// author: Rhys Whitley
// institution: Macquarie University
// date: Jan 2015

class Energy {
/* Holds all calculations for the transfer of energy between the atmosphere and
 * the earth's surface and the describes the conversion of energy to water
 * vapour
 */
    private:
        // numerical constants used in evaporation calculations
        static const float
                    eccn    = 0.01675,
                    d2r     = M_PI/180.0,
                    cw      = 1.0,
                    solc    = 1367.0,
                    dtime   = 0.0036,
                    nsecs   = 86400.0,
                    albedo  = 0.17;

    public:
        // single returns
        float slope_vapor_temperature   ( const float tair  );
        float solar_const               ( const int day     );
        float change_solar_angle        ( const int day     );
        float terrestrial_radiation     ( const int day,    const float lat );
        float par_atmos_top             ( const float msolc,const float hs, const float x,      const float y );
        float actual_evaporation        ( const float par,  const float h1, const float spl,    const float u,  const float v,  const float dcon );
        float equilibrium_evaporaiton   ( const float h0,   const float u,  const float v );
        float photo_active_radiation    ( const float rabs, const float hs, const float rsin,   const float rcos );
        float calc_parameter_hs         ( const float x,    const float y,  const float lat,    const float delta   );
        // multiple returns
        void calc_parameters_h01        ( const float u,    const float v,  const float spl,    float &h0,      float &h1 );
        void euclidean_projection       ( const float lat,  const float delta,  float& x,       float& y);
        void energy_lookup_table        ( const float temp, float &gamma,       float &lambda );
        void sunshine                   ( const float swr,  const float lat,    float& sunhr,   float& totrad );
        // main routines
        // point
        void evaporate                  ( const float tair, const float fsun,   const float lat, const float spl, const int day,
                                                float& par, float& aet, float& eet, float& pet, float& det );
        // gridded
        void grid_evaporate_wrapper     ( GridCell &gc,     const float spl,    const int day );
        void waterBucket                ( GridCell &gc );
        // inline functions
        inline float solar_angle        ( const int day )                                   { return 360.*day/365; }
        inline float sunlight_hours     ( const float swr, const float rtoa )               { return 2*(swr/rtoa-0.25); }
        inline float absorb_radiation   ( const float fsun )                                { return (0.25+0.5*fsun)*(1-albedo); }
        inline float reflect_radiation  ( const float fsun,const float tair )               { return (0.2+0.8*fsun)*(107.-tair); }
        inline float energy_convert     ( const float lam, const float gam, const float dS ){ return dtime/lam*dS/(dS+gam); }
    };
#endif
