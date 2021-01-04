#include <iostream>
#include <cmath>
#include <string>
#include "spectrum2rgb.h"
using namespace std;

// wavelengths are always in nm; or in angstrom and then converted to nm.
int main(int argc, char* argv[]){
  FILE *fpin;
  FILE *fpout;		
				
  fpin = stdin;
  fpout = stdout;
	
	float multiplication_factor = 1.0;
		
  float wavelength;
	float intensity;
	
	bool BB_flag = false;
	float BB_T;
	double BB_intensity;
	float BB_wavelength;
	
	bool clip_flag = false;
	
	float response_functions_values[] = {0.0, 0.0, 0.0};
	
	double X = 0.0;
	double Y = 0.0;
	double Z = 0.0;
	
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	
	double norm;
	
	double r = 0.0;
	double g = 0.0;
	double b = 0.0;
	
	
	float normalisationfactor = 1.0;
	
/*	
	cout << S_inv[0][0] << " " << S_inv[0][1] << " " << S_inv[0][2] << endl;
	cout << S_inv[1][0] << " " << S_inv[1][1] << " " << S_inv[1][2] << endl;
	cout << S_inv[2][0] << " " << S_inv[2][1] << " " << S_inv[2][2] << endl <<
	endl;
	
S_inv[0][0] = S_inv_srgb[0][0];
S_inv[0][1] = S_inv_srgb[0][1];
S_inv[0][2] = S_inv_srgb[0][2];

S_inv[1][0] = S_inv_srgb[1][0];
S_inv[1][1] = S_inv_srgb[1][1];
S_inv[1][2] = S_inv_srgb[1][2];

S_inv[2][0] = S_inv_srgb[2][0];
S_inv[2][1] = S_inv_srgb[2][1];
S_inv[2][2] = S_inv_srgb[2][2];

	cout << S_inv[0][0] << " " << S_inv[0][1] << " " << S_inv[0][2] << endl;
	cout << S_inv[1][0] << " " << S_inv[1][1] << " " << S_inv[1][2] << endl;
	cout << S_inv[2][0] << " " << S_inv[2][1] << " " << S_inv[2][2] << endl <<
	endl;
*/

	for(int i = 0; i < argc; i++) {
			if (!strcmp(argv[i],"-a")) {
				multiplication_factor = 0.1;	
			}				
			if (!strcmp(argv[i],"-b")) {
					BB_flag = true;
					BB_T = atof(argv[i + 1]); 
			}		
			if (!strcmp(argv[i],"-c")) {
					clip_flag = true;
			}					
			if (!strcmp(argv[i],"-f")) {
					fpin = fopen (argv[i + 1], "rt"); 
			}
			if (!strcmp(argv[i],"-o")) {
					fpout = fopen (argv[i + 1], "wt"); 
			}			
			if (!strcmp(argv[i],"-w")) {
					xw =  atof(argv[i + 1]);
					yw =  atof(argv[i + 2]);
			}			
			if (!strcmp(argv[i],"-pr")) {
					xr =  atof(argv[i + 1]);
					yr =  atof(argv[i + 2]);
			}		
			if (!strcmp(argv[i],"-pg")) {
					xg =  atof(argv[i + 1]);
					yg =  atof(argv[i + 2]);
			}								
			if (!strcmp(argv[i],"-pb")) {
					xb =  atof(argv[i + 1]);
					yb =  atof(argv[i + 2]);
			}		
			if (!strcmp(argv[i],"-help") || !strcmp(argv[i],"--help") || !strcmp(argv[i],"-h")) {
					cerr << "Calculate the RGB components of an optical spectrum or a black   " << endl;
					cerr << "body of a certain temperature. The spectrum can be fed as a two  " << endl;
					cerr << "dimensional list of intensity vs. wavelength.                    " << endl << endl;
					cerr << "Options: -a          spectrum in angstrom instead of nm (default)" << endl;
					cerr << "         -b T        give the rgb values for a black body with   " << endl;
					cerr << "                     temperature T  0.03 K < T < 10^20 K            " << endl;
					cerr << "         -c clip     the rgb values at 0 and 1                   " << endl;
 					cerr << "         -f file     send output to file                         " << endl;
					cerr << "         -o file     get the spectrum from file                  " << endl;
					cerr << "         -w xw yw    adjust the whitepoint                       " << endl;
					cerr << "         -pr xr yr   adjust the red primary                      " << endl;
					cerr << "         -pg xg yg   adjust the green primary                    " << endl;
					cerr << "         -pb xb yb   adjust the blue primary                     " << endl << endl;	
					return 0;														
			}		
	}
		

// some constants, useful later		
float Ca = ((yg - yb)*(xw - xb) + (xb - xg)*(yw - yb))/((xr - xb)*(yg - yb) - (xg - xb)*(yr - yb));
float Cb = ((yb - yr)*(xw - xb) + (xr - xb)*(yw - yb))/((xr - xb)*(yg - yb) - (xg - xb)*(yr -
yb));
float Cc = 1.0 - Ca - Cb;

// determine the conversionmatrix' components
float S_inv[3][3] = {
{((xg - 1)*yb + yg - xb*yg)/(Ca*(xr*(yg -yb) + xg*(yb - yr) + xb*(yr - yg))),

 (xg - xg*yb + xb*(yg -1))/(Ca*(xr*(yb - yg) + xb*(yg - yr) + xg*(yr - yb))),
 
 (xg*yb - xb*yg)/(Ca*(xr*(yg - yb) + xg*(yb - yr) + xb*(yr - yg)))},
 
 
{((xr - 1)*yb + yr - xb*yr)/(Cb*(xr*(yb - yg) + xb*(yg - yr) + xg*(yr - yb))),

 (xr - xr*yb + xb*(yr - 1))/(Cb*(xr*(yg - yb) + xg*(yb - yr) + xb*(yr - yg))),
 
 (xr*yb - xb*yr)/(Cb*(xr*(yb - yg) + xb*(yg - yr) + xg*(yr - yb)))},
 
 
{((xr - 1)*yg + yr - xg*yr)/(Cc*(xr*(yg - yb) + xg*(yb - yr) + xb*(yr - yg))),

 (xr - xr*yg + xg*(yr - 1))/(Cc*(xr*(yb - yg) + xb*(yg - yr) + xg*(yr - yb))),
 
 (xr*yg - xg*yr)/(Cc*(xr*(yg - yb) + xg*(yb - yr) + xb*(yr - yg)))}
};
		
	
	if (BB_flag) {
	//use a (normalised) Black Body curve
	// 25 < T < 10^20
		
		for (int i = 0; i < CIE_1964_response_functions_length; i++) {
			
			BB_wavelength = CIE_1964_response_functions_min + i;
			
			BB_intensity = get_BB_norm_function_value(BB_wavelength, BB_T);

			get_response_function_values(BB_wavelength, response_functions_values);
			
			//sum the response functions multiplied with the spectrum to obtain X, Y, Z
			
			X += response_functions_values[0]*BB_intensity;
			Y += response_functions_values[1]*BB_intensity;			
			Z += response_functions_values[2]*BB_intensity;
			
		}
	}
	else {	
		//use a spectrum from file or stdin, or where-ever
    while( !feof(fpin) ){
			fscanf(fpin,"%f %f", &wavelength, &intensity);
			
			wavelength*=multiplication_factor;
			
			get_response_function_values(wavelength, response_functions_values);
			
			//sum the response functions multiplied with the spectrum to obtain X, Y, Z
			
			X += response_functions_values[0]*intensity;
			Y += response_functions_values[1]*intensity;			
			Z += response_functions_values[2]*intensity;
		}
	}
		
	//the X, Y, Z values are now known
	//normalise these to get x, y, z
	norm = 1.0/(X + Y + Z);
	x = X*norm;
	y = Y*norm;
	z = Z*norm;
		
	//let matrix S_inv work on x, ,y, z;
		
	r = S_inv[0][0]*x + S_inv[0][1]*y + S_inv[0][2]*z;
	g = S_inv[1][0]*x + S_inv[1][1]*y + S_inv[1][2]*z	;			
	b = S_inv[2][0]*x + S_inv[2][1]*y + S_inv[2][2]*z;

	if (clip_flag){
		//clip the colours which lie outside the gamut
		if (r > 1.0)
			r = 1.0;	
		if (r < 0.0)
			r = 0.0;
		if (g > 1)
			g = 1.0;
		if (g < 0.0)
			g = 0.0;
		if (b > 1)
			b = 1.0;
		if (b < 0.0)
			b = 0.0;
	}
	else{
		// alternative adjustment where the rgb values are all made positive
		// and (later) normalised
		float add_value = 0.0;
		
		if (r < 0.0)
			add_value += 0.0 - r;
		if (g < 0.0)
			add_value += 0.0 - g;
		if (b < 0.0)
			add_value += 0.0 - b;
		
		r += add_value;
		g += add_value;
		b += add_value;		
	}
	
	//normalise the rgb values such that one component equals unity
	normalisationfactor = r;
	if (g > r) {
		normalisationfactor = g;
	}
	if (b > normalisationfactor) {
			normalisationfactor = b;
	}
	
	normalisationfactor = 1.0/normalisationfactor;
	
	r *= normalisationfactor;
	g *= normalisationfactor;
	b *= normalisationfactor;	
	

 	fprintf(fpout, "%f %f %f\n",r, g, b);	
		
	return 0;
}
/*
if the rgb values obtained are > 1 or < 0, then they lie
outside the chosen colour gamut.
*/

void get_response_function_values(float wavelength, float* return_array) {
	
	if ( (wavelength >= CIE_1964_response_functions[0][0]) &&
			 (wavelength <= CIE_1964_response_functions[CIE_1964_response_functions_length][0]) ) {
		
		float sample_width = (CIE_1964_response_functions_max -
				 CIE_1964_response_functions_min)/(CIE_1964_response_functions_length);
		
		int closest_lower_index = (int)((wavelength-CIE_1964_response_functions_min)
				                                                           /sample_width);
		
		float lower_wavelength = CIE_1964_response_functions[closest_lower_index][0];
		float upper_wavelength = CIE_1964_response_functions[closest_lower_index + 1][0];
		
		float relative_difference = wavelength - lower_wavelength;
		
		float interpolation_factor = relative_difference/(upper_wavelength - lower_wavelength);
				
				
		// the three response functions values
		
		return_array[0] = CIE_1964_response_functions[closest_lower_index][1] + 
				(CIE_1964_response_functions[closest_lower_index + 1][1] -
				 CIE_1964_response_functions[closest_lower_index][1])*interpolation_factor;
				
		return_array[1] = CIE_1964_response_functions[closest_lower_index][2] + 
				(CIE_1964_response_functions[closest_lower_index + 1][2] -
				 CIE_1964_response_functions[closest_lower_index][2])*interpolation_factor;
				 				
		return_array[2] = CIE_1964_response_functions[closest_lower_index][3] + 
				(CIE_1964_response_functions[closest_lower_index + 1][3] -
				 CIE_1964_response_functions[closest_lower_index][3])*interpolation_factor;
		
	}
	else {
		return_array[0] = 0.0;
		return_array[1] = 0.0;
		return_array[2] = 0.0;		
	}
	
}

double get_BB_norm_function_value(float wavelength, float T) {
	//wavelength in nm
	//constants in SI
	
	double lambda = wavelength*pow(10.0,-9.0); //lambda_var is in nm, so convert to m
	double hcoverk = 0.01437605113;
	double  lambda_naught = 0.002900237/T;
	double B = pow((lambda_naught/lambda),5)*
			(exp(hcoverk/(lambda_naught*T)) - 1)/(exp(hcoverk/(lambda*T)) - 1);
	
	/*
		Here numerical precision of a float will fail.
		Use the long wavlength approximation and scale
		the exponent to the longest wavelength we have.
	*/
	if (T < 25) {
		B = pow((lambda_naught/lambda),5)*
			(exp(hcoverk/(lambda_naught*T)) - 1)*exp(-hcoverk/(lambda*T) +
				 hcoverk/(CIE_1964_response_functions_max*pow(10.0,-9.0)*T));
	}			
	
	return B;
}
