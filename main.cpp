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
	// speed of light in microns per second
	Doub C = 2.998e14;
	// number of bands
	Int N_band = 5;

	// These are used for allocating memory, could cause segfaults maybe
	// maximum number of initial masses being read in from IMF
	Int maxN_IMF = 6000000;
	// maximum number of lines to be read from isochrone file
	Int maxN_iso = 1600;

	// All below for Z, Y, J, H & F
	// Conversion from AB to Vega magnitudes 
	Doub AB_Vega[] = {0.487, 0.653, 0.958, 1.287, 1.552}; // This is m_AB - m_Vega
	// AB zeropoints from arXiv:1702.01747
	Doub ab_zeros[] = {26.39, 26.41, 26.35, 26.41, 25.96};
	// central wavelength in microns
	Doub lam_c[] = {0.87,1.09,1.30,1.60,1.88};
	// central frequency in Hz
	Doub nu_c[] = {0.,0.,0.,0.,0.};
	Int j;
	for(j=0;j<N_band;j++){
		nu_c[j] = C/lam_c[j];
	}

	/* 5 sigma point source detection limiting magnitudes for a 1000 second exposure
	* in each band, calculated from Chris Hirata's ETC code, ranges for bands used are
	* R = 0.5-0.8 microns,   Z = 0.75-1.0 microns,   Y = 0.9-1.2 microns
	* J = 1.1-1.5 microns,   W = 0.95-2.0 microns,   H = 1.35-1.8 microns, F = 1.65-2.0 microns 
	* Z, Y, J, H, F */
	Doub ps_detect_5slim[] = {26.061,25.989,25.939,25.819,25.517};
	for(j=0;j<N_band;j++){
		ps_detect_5slim[j] += 2.5 *log10(sqrt(t_expose/1000.));
	}

	// AB zero-point flux definitin in ergs s^-1 Hz^-1
	Doub f_ab_zero = 6.626e-27;


	/* Read in initializing files such as the IMF 
	 * and the isochrones 
	 */
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

	/* Reading in the isochrone so that the IMF can be 
	 * interpolated and turned in to colors/magnitudes
	 */
	Doub *imass, *Z, *Y, *J, *H, *F;
	imass = (Doub*) malloc(sizeof(Doub)*maxN_iso);
	Z = (Doub*) malloc(sizeof(Doub)*maxN_iso);
	Y = (Doub*) malloc(sizeof(Doub)*maxN_iso);
	J = (Doub*) malloc(sizeof(Doub)*maxN_iso);
	H = (Doub*) malloc(sizeof(Doub)*maxN_iso);
	F = (Doub*) malloc(sizeof(Doub)*maxN_iso);

	// Reading isochrone file
	ifstream isofile ("samp_iso.txt");
	if(!isofile.is_open()){
		cerr << "There was a problem opening the isochrone file \n";
		return 1;
	}
	Int i = 0;
	Doub min_imass=10., max_imass=0.;
	while (isofile>>temp) {
		imass[i] = temp;
		if (imass[i]>max_imass)
			max_imass = imass[i];
		if (imass[i]<min_imass)
			min_imass = imass[i];
		isofile>>temp;
		Z[i] = temp;
		isofile>>temp;
		Y[i] = temp;
		isofile>>temp;
		J[i] = temp;
		isofile>>temp;
		H[i] = temp;
		isofile>>temp;
		F[i] = temp;
		i += 1;
	}
	isofile.close();



	// de-allocate memory
	free(kmass);
	free(imass);
	free(Z);
	free(Y);
	free(J);
	free(H);
	free(F);

	/* code */
	return 0;
}
