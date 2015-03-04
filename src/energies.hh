#ifndef _energies_H
#define _energies_H

// author: Rhys Whitley
// institution: Macquarie University
// date: Mar 2015

class Energies {
/* Holds all calculations for the transfer of energy between the atmosphere and
 * the earth's surface and the describes the conversion of energy to water
 * vapour
 */
    private:

        float   aet,    // actual evaporation
                pet,    // potential evaporation
                eet,    // equilibrium evaporation
                det,    // difference between actual and equilibrium evaporation
                ppt,    // rainfall
                run,    // runoff or overland flow
                swc,    // water content of the soil profile
                par,    // photosynthetically active radiation
                fsun,   // fractional sunshine hours
                srad,   // incoming solar radiation
                rtoa,   // radiation at the top of the atmosphere
                tmin,   // minimum ambient air temperature
                tmax,   // maximum ambient air temperature
                tmean,  // mean ambient air temperature
                trange; // range of ambient air temperature

    public:
        // getters and setters for manipulating energy variables
        void set_aet    ( const float value );
        void set_pet    ( const float value );
        void set_eet    ( const float value );
        void set_det    ( const float value );
        void set_ppt    ( const float value );
        void set_run    ( const float value );
        void set_swc    ( const float value );
        void set_par    ( const float value );
        void set_fsun   ( const float value );
        void set_srad   ( const float value );
        void set_rtoa   ( const float value );
        void set_tmin   ( const float value );
        void set_tmax   ( const float value );
        void set_tmean  ( const float value );
        void set_trange ( const float value );

        float get_aet    () const;
        float get_pet    () const;
        float get_eet    () const;
        float get_det    () const;
        float get_ppt    () const;
        float get_run    () const;
        float get_swc    () const;
        float get_par    () const;
        float get_fsun   () const;
        float get_srad   () const;
        float get_rtoa   () const;
        float get_tmin   () const;
        float get_tmax   () const;
        float get_tmean  () const;
        float get_trange () const;

        // declaration of other functions
        void reset_all_values( const float value );

    };
#endif

