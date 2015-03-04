#ifndef _tempVar_H
#define _tempVar_H

// author: Rhys Whitley
// institution: Macquarie University
// date: Mar 2015

class ClimateIndex {
/* Holds all calculations for the temperature and sensible heat
 */
    private:

        unsigned
        int     gdd0,   // growing degree days: number of days above 0 Celsius
                gdd5,   // growing degree days: number of days above 5 Celsius
                gdd10,  // growing degree days: number of days above 10 Celsius
                chill;  // chilly days: number of days below 0 Celsius
        float   alpha,  // ratio of actual evaporation to equilibrium evaporation
                moist;  // ratio of precipitation to potential evaporation

    public:
        // getters and setters for manipulating energy variables
        void set_gdd0   ( const float value );
        void set_gdd5   ( const float value );
        void set_gdd10  ( const float value );
        void set_chill  ( const float value );
        void set_alpha  ( const float value );
        void set_moist  ( const float value );

        float get_gdd0  () const;
        float get_gdd5  () const;
        float get_gdd10 () const;
        float get_chill () const;
        float get_alpha () const;
        float get_moist () const;

        // declaration of other functions
        void reset_all_values( const float value );

    };
#endif

