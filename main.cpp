/* Code for calculating the expected number of resolved stars 
 * for a given surface brightness (mu), distance to the observed
 * object (d), stellar population, and exposure time (t)
 * 
 * This implementation is in C++ just to make sure that this 
 * calculation is quick and able to be used easily by the community.
 *
 * written by Lachlan Lancaster, March 2020
 */

#include <fstream>
#include <cmath>
#include <iostream>
#include <array>
#include <algorithm>
#include <random>
#include "utils.h"

using namespace std;


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

	// defintion of field of view
	double sqdeg_sqas = 3600.*3600.;
	double dra=1.,ddec = 1.;
	double FoV_as = dra*ddec*sqdeg_sqas;


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
	double expose_fac = 2.5*log10(sqrt(t_expose/1000.));
	for(j=0;j<N_band;j++){
		ps_detect_5slim[j] += expose_fac;
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

	// interpolate on to the IMF and convert from
	// Vega magnitudes to AB magnitudes
	int N_imf_filt=0;
	for(j=0;j<N_IMF;j++){
		if ((kmass[j] < max_imass) && (kmass[j]>min_imass)){
			Z_out[N_imf_filt] = linear_interp(kmass[j],imass,Z) + AB_Vega[0] + mu;
			Y_out[N_imf_filt] = linear_interp(kmass[j],imass,Y) + AB_Vega[1] + mu;
			J_out[N_imf_filt] = linear_interp(kmass[j],imass,J) + AB_Vega[2] + mu;
			H_out[N_imf_filt] = linear_interp(kmass[j],imass,H) + AB_Vega[3] + mu;
			F_out[N_imf_filt] = linear_interp(kmass[j],imass,F) + AB_Vega[4] + mu;
			N_imf_filt++;
		}
	}

	// allocate a index array that will be repeatedley shuffled
	// in order to Monte-Carlo sample the IMF
	array<int,10> index_arr = {0,1,2,3,4,5,6,7,8,9};
	/*for (i=0;i<4;i++){
		index_arr[i] = i;
		}*/
	unsigned seed = 42;
	shuffle (index_arr.begin(), index_arr.end(), default_random_engine(seed));

	//cout << "Shuffled elements: ";
	//for(int &x: index_arr) cout << ' ' << x;
	//cout << "\n";


	
	// set uppper flux limits
	double S_J = S + AB_Vega[2] - ab_zeros[2];
	double f_S_J = f_ab_zero*nu_c[2] * pow(10,((S_J)/(-2.5)));
	double f_S_J_tot = f_S_J*FoV_as;
	double S_H = S + AB_Vega[3] - ab_zeros[3];
	double f_S_H = f_ab_zero*nu_c[3] * pow(10,((S_H)/(-2.5)));
	double f_S_H_tot = f_S_H*FoV_as;
	double S_F = S + AB_Vega[4] - ab_zeros[4];
	double f_S_F = f_ab_zero*nu_c[4] * pow(10,((S_F)/(-2.5)));
	double f_S_F_tot = f_S_F*FoV_as;

	cout << "Flux limit in J band is " << f_S_J_tot << "ergs/s" << "\n";

	double *f_Z, *f_Y, *f_J, *f_H, *f_F;
	f_Z = (double *) malloc(sizeof(double)*N_imf_filt);
	f_Y = (double *) malloc(sizeof(double)*N_imf_filt);
	f_J = (double *) malloc(sizeof(double)*N_imf_filt);
	f_H = (double *) malloc(sizeof(double)*N_imf_filt);
	f_F = (double *) malloc(sizeof(double)*N_imf_filt);


	//fluxes in photon counts per second
	for (i=0;i<N_imf_filt;i++){
		f_Z[i] = f_ab_zero*nu_c[0] * pow(10, (Z_out[i] - ab_zeros[0])/(-2.5));
		f_Y[i] = f_ab_zero*nu_c[1] * pow(10, (Y_out[i] - ab_zeros[1])/(-2.5));
		f_J[i] = f_ab_zero*nu_c[2] * pow(10, (J_out[i] - ab_zeros[2])/(-2.5));
		f_H[i] = f_ab_zero*nu_c[3] * pow(10, (H_out[i] - ab_zeros[3])/(-2.5));
		f_F[i] = f_ab_zero*nu_c[4] * pow(10, (F_out[i] - ab_zeros[4])/(-2.5));
	}


	// cumulative sum of flux array
	double CS_fJ = 0;
	int n_needed = 0;
	while((CS_fJ < f_S_J_tot) && (n_needed < N_imf_filt)){
		CS_fJ += f_J[n_needed];
		n_needed++;
	}
	// check to make cure that the whole filter wasn't needed
	if (n_needed==N_imf_filt){
		cout << "Used the entire IMF. You should decrease the size of the field or increase the IMF sample." << "\n";
	}

	// count the number that are above the detection threshold
	int n_detected = 0;
	for(i=0;i<n_needed;i++){
		if(J_out[i]<ps_detect_5slim[2]){
			n_detected++;
		}
	}
	cout << n_needed << " stars are in the field of view" << "\n";
	cout << n_detected << " stars are detected per square degree" << "\n";


	// de-allocate memory
	free(kmass);
	free(imass);
	free(Z); free(Z_out); free(f_Z);
	free(Y); free(Y_out); free(f_Y);
	free(J); free(J_out); free(f_J);
	free(H); free(H_out); free(f_H);
	free(F); free(F_out); free(f_F);

	/* code */
	return 0;
}
