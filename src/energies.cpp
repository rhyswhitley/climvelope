#include "energies.hh"

void Energies::reset_all_values( const float value ) {
// sets all energy variables to some default value
    par     = value;
    aet     = value;
    eet     = value;
    pet     = value;
    det     = value;
    ppt     = value;
    run     = value;
    fsun    = value;
    srad    = value;
    rtoa    = value;
    tmin    = value;
    tmax    = value;
    tmean   = value;
    trange  = value;
}

// Setters
void Energies::set_par     ( const float value ) { par = value; }
void Energies::set_aet     ( const float value ) { aet = value; }
void Energies::set_pet     ( const float value ) { pet = value; }
void Energies::set_eet     ( const float value ) { eet = value; }
void Energies::set_det     ( const float value ) { det = value; }
void Energies::set_ppt     ( const float value ) { ppt = value; }
void Energies::set_run     ( const float value ) { run = value; }
void Energies::set_swc     ( const float value ) { swc = value; }
void Energies::set_fsun    ( const float value ) { fsun = value; }
void Energies::set_srad    ( const float value ) { srad = value; }
void Energies::set_rtoa    ( const float value ) { rtoa = value; }
void Energies::set_tmin    ( const float value ) { tmin = value; }
void Energies::set_tmax    ( const float value ) { tmax = value; }
void Energies::set_tmean   ( const float value ) { tmean = value; }
void Energies::set_trange  ( const float value ) { trange = value; }

// Getters
float Energies::get_par    () const { return par; }
float Energies::get_aet    () const { return aet; }
float Energies::get_pet    () const { return pet; }
float Energies::get_eet    () const { return eet; }
float Energies::get_det    () const { return det; }
float Energies::get_ppt    () const { return ppt; }
float Energies::get_run    () const { return run; }
float Energies::get_swc    () const { return swc; }
float Energies::get_fsun   () const { return fsun; }
float Energies::get_srad   () const { return srad; }
float Energies::get_rtoa   () const { return rtoa; }
float Energies::get_tmin   () const { return tmin; }
float Energies::get_tmax   () const { return tmax; }
float Energies::get_tmean  () const { return tmean; }
float Energies::get_trange () const { return trange; }

