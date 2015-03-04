#include "energies.hh"
#include "climindex.hh"
#include <iostream>
#include <vector>

using namespace std;

int main() {

    vector<Energies> eVars(10);
    vector<ClimateIndex> eInd(10);

    for( unsigned int i=0; i<eVars.size(); i++ ){

        eVars[i].reset_all_values(-999.);
        eVars[i].set_aet(1000.);

        eInd[i].reset_all_values(-999.);
        eInd[i].set_gdd0(100);

        cout << eVars[i].get_aet() << endl;
        cout << eInd[i].get_gdd0() << endl;
    }

    cout << "Hey it worked!" << endl;

    return 0;
}
