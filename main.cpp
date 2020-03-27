/* Code for calculating the expected number of resolved stars 
 * for a given surface brightness (mu), distance to the observed
 * object (d), stellar population, and exposure time (t)
 * 
 * This implementation is in C++ just to make sure that this 
 * calculation is quick and able to be used easily by the community.
 *
 * written by Lachlan Lancaster, March 2020
 */

#include <omp.h>
#include "nr3.h"
#include "utils.h"

int main()
{
	/* definition of global parameters, probably will 
	 * be read in from stdin in later development 
	 */

	// surface birghtness in mags/arcsec^2
	Doub S = 33.; 
	// exposure time in seconds
	Doub t_expose = 1000.;
	// distance in Mpc
	Doub d = 1.;
	// distance modulus 
	Doub mu = 5.*log10(d*1e5);

	/* Read in initializing files such as the IMF 
	 * and the isochrones 
	 */
	// maximum number of initial masses being read in
	Int maxN_IMF = 6000000;
	Doub temp;
	Doub *kmass;
	kmass = (Doub*) malloc(sizeof(Doub)*maxN_IMF);

	// reading in IMF file
	ifstream massfile ("mid_set.txt");
	if(!massfile.is_open()){
		cerr << "There was a problem opening the IMF file \n";
		return 1;
	}
	Int i = 0;
	while (massfile>>temp) {
		kmass[i] = temp;
		i += 1;
	}
	massfile.close();

	// 


	/* code */
	return 0;
}