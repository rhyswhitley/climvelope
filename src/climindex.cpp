#include "climindex.hh"

void ClimateIndex::reset_all_values( const float value ) {
// sets all energy variables to some default value
    gdd0    = int(value);
    gdd5    = int(value);
    gdd10   = int(value);
    chill   = value;
    alpha   = value;
    moist   = value;
}

// Setters
void ClimateIndex::set_gdd0  ( const float value ) { gdd0 = value; }
void ClimateIndex::set_gdd5  ( const float value ) { gdd5 = value; }
void ClimateIndex::set_gdd10 ( const float value ) { gdd10 = value; }
void ClimateIndex::set_chill ( const float value ) { chill = value; }
void ClimateIndex::set_alpha  ( const float value ) { alpha = value; }
void ClimateIndex::set_moist  ( const float value ) { moist = value; }

// Getters
float ClimateIndex::get_gdd0     () const { return gdd0; }
float ClimateIndex::get_gdd5     () const { return gdd5; }
float ClimateIndex::get_gdd10    () const { return gdd10; }
float ClimateIndex::get_chill    () const { return chill; }
float ClimateIndex::get_alpha  () const { return alpha; }
float ClimateIndex::get_moist  () const { return moist; }

