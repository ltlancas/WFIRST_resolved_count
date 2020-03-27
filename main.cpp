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
	double S = 33.; 
	// exposure time in seconds
	double t_expose = 1000.;
	// distance in Mpc
	double d = 1.;
	// distance modulus 
	double mu = 5.*log10(d*1e5);
	// speed of light in microns per second
	double C = 2.998e14;
	// number of bands
	int N_band = 5;

	// These are used for allocating memory, could cause segfaults maybe
	// maximum number of initial masses being read in from IMF
	int maxN_IMF = 6000000;
	// maximum number of lines to be read from isochrone file
	int maxN_iso = 1600;

	// All below for Z, Y, J, H & F
	// Conversion from AB to Vega magnitudes 
	double AB_Vega[] = {0.487, 0.653, 0.958, 1.287, 1.552}; // This is m_AB - m_Vega
	// AB zeropoints from arXiv:1702.01747
	double ab_zeros[] = {26.39, 26.41, 26.35, 26.41, 25.96};
	// central wavelength in microns
	double lam_c[] = {0.87,1.09,1.30,1.60,1.88};
	// central frequency in Hz
	double nu_c[] = {0.,0.,0.,0.,0.};
	int j;
	for(j=0;j<N_band;j++){
		nu_c[j] = C/lam_c[j];
	}

	/* 5 sigma point source detection limiting magnitudes for a 1000 second exposure
	* in each band, calculated from Chris Hirata's ETC code, ranges for bands used are
	* R = 0.5-0.8 microns,   Z = 0.75-1.0 microns,   Y = 0.9-1.2 microns
	* J = 1.1-1.5 microns,   W = 0.95-2.0 microns,   H = 1.35-1.8 microns, F = 1.65-2.0 microns 
	* Z, Y, J, H, F */
	double ps_detect_5slim[] = {26.061,25.989,25.939,25.819,25.517};
	for(j=0;j<N_band;j++){
		ps_detect_5slim[j] += 2.5 *log10(sqrt(t_expose/1000.));
	}

	// AB zero-point flux definitin in ergs s^-1 Hz^-1
	double f_ab_zero = 6.626e-27;


	/* Read in initializing files such as the IMF 
	 * and the isochrones 
	 */
	int N_IMF;
	double temp;
	double *kmass;
	kmass = (double*) malloc(sizeof(double)*maxN_IMF);

	// reading in IMF file
	ifstream massfile ("mid_set.txt");
	if(!massfile.is_open()){
		cerr << "There was a problem opening the IMF file \n";
		return 1;
	}
	int i = 0;
	while (massfile>>temp) {
		kmass[i] = temp;
		i += 1;
	}
	N_IMF = i;
	massfile.close();

	/* Reading in the isochrone so that the IMF can be 
	 * interpolated and turned in to colors/magnitudes
	 */
	int N_iso;
	double *imass, *Z, *Y, *J, *H, *F;
	double *Z_out, *Y_out, *J_out, *H_out, *F_out;
	imass = (double*) malloc(sizeof(double)*maxN_iso);
	Z = (double*) malloc(sizeof(double)*maxN_iso);
	Y = (double*) malloc(sizeof(double)*maxN_iso);
	J = (double*) malloc(sizeof(double)*maxN_iso);
	H = (double*) malloc(sizeof(double)*maxN_iso);
	F = (double*) malloc(sizeof(double)*maxN_iso);
	Z_out = (double*) malloc(sizeof(double)*N_IMF);
	Y_out = (double*) malloc(sizeof(double)*N_IMF);
	J_out = (double*) malloc(sizeof(double)*N_IMF);
	H_out = (double*) malloc(sizeof(double)*N_IMF);
	F_out = (double*) malloc(sizeof(double)*N_IMF);

	// Reading isochrone file
	ifstream isofile ("samp_iso.txt");
	if(!isofile.is_open()){
		cerr << "There was a problem opening the isochrone file \n";
		return 1;
	}
	i = 0;
	double min_imass=10., max_imass=0.;
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
	N_iso = i;
	isofile.close();

	int N_imf_filt=0;
	for(j=0;j<N_IMF;j++){
		if ((kmass[j] < max_imass) && (kmass[j]>min_imass)){
			Z_out[N_imf_filt] = linear_interp(kmass[j],imass,Z);
			Y_out[N_imf_filt] = linear_interp(kmass[j],imass,Y);
			J_out[N_imf_filt] = linear_interp(kmass[j],imass,J);
			H_out[N_imf_filt] = linear_interp(kmass[j],imass,H);
			F_out[N_imf_filt] = linear_interp(kmass[j],imass,F);
			N_imf_filt++;
		}
	}


	// de-allocate memory
	free(kmass);
	free(imass);
	free(Z); free(Z_out);
	free(Y); free(Y_out);
	free(J); free(J_out);
	free(H); free(H_out);
	free(F); free(F_out);

	/* code */
	return 0;
}
