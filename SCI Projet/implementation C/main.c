//====================================================================//
// Code property:                                                     //
// ==============                                                     //
// TELECOM BRETAGNE                                                   //
// Dpt. Signal et Communications                                      //
// Technopole Brest-Iroise                                            //
// CS 83818 - 29238 Brest Cedex 3 - France                            //
//                                                                    //
// General info:                                                      //
// =============                                                      //
// program: SOBI.c                                                    //
// Adapted from: program SeparationSourceSobi.c, 26/09/2013,          //
// by J. Trubuil (SC)                                                 //
// last upate: SOBI.c, 02/11/2016                                     //
// info: thierry.chonavel@telecom-bretagne.eu                         //
//       thiery.legall1@telecom-bretagne.eu                           //
//                                                                    //
// Code info:  	                                                      //
// ==========                                                         //
// Separation of a two signals mixture observed on two sensors        //
// using SOBI algorithm                                               //
// Mixing matrix: A=[1,1;-1,2]                                        //
// Mixed signals files:  sobi_in_1.txt & sobi_in_2.txt                //
// (generated with Matlab program  "CreationMelange.m")               //
// Outputs: files sobi_out_1.txt & sobi_out_2.txt                     //
//                                                                    //
// Source code compilation and execution:                             //
// ======================================                             //
// under Ubuntu, compile and link with libm.so as follows             //
// note : replace separation.* by separation_v0.* for initial test    //
//                                                                    //
// $gcc -c math_functions.c -o math_functions.o                       //
// $gcc -c separation.c -o separation.o                               //
// $gcc SOBI.c -lm math_functions.o separation.o -o ./SOBI.exe        //
// $chmod +x SOBI.exe                                                 //
// $./SOBI.exe                                                        //
//====================================================================//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include "math_functions.h"
#include "separation_v0.h"
//#include "separation.h"

int main ()
{
  // Declaration of file pointers
  FILE * file_input_1,     * file_input_2;
  FILE * file_output_1,    * file_output_2;
  FILE * file_reference_1, * file_reference_2;

  // Declaration of array pointers
  double * input_1      = NULL;
  double * input_2      = NULL;
  double * output_1     = NULL;
  double * output_2     = NULL;
  double * reference_1  = NULL;
  double * reference_2  = NULL;

  // Declaration of variables
  int    nb_samples, n1, n2;          // number of samples
  int    time_delay     = 10;          // number of time delays for time correlations in SOBI
  int    i ,j;                        // loop variables
  double cos_1, cos_2, err_1, err_2;  // performance indices

  printf("*******************************************\n");
  printf("* Sources separation with SOBI algorithm  *\n");
  printf("*******************************************\n");
	// Openning input files
	file_input_1 = fopen("sobi_in_1.txt","rb");
        file_input_2 = fopen("sobi_in_2.txt","rb");
	if(file_input_1==NULL || file_input_2==NULL) {  // Control files openning
		printf("Erreur openning files for reading");
		return EXIT_FAILURE;
	}
	// Openning output files
  	file_output_1 = fopen("sobi_out_1.txt","wb");
  	file_output_2 = fopen("sobi_out_2.txt","wb");
	if(file_output_1==NULL || file_output_2==NULL) { // Control files openning
		printf("Erreur openning files for writting");
    		return EXIT_FAILURE;
  	}

	// Number of samples to be processed
	fseek(file_input_1,0L,SEEK_END);         // position indices at last files position
  	fseek(file_input_2,0L,SEEK_END);
	n1 = ftell(file_input_1)/sizeof(double); // compute the number of doubles in the files
	n2 = ftell(file_input_1)/sizeof(double);
	if(n1 < n2) nb_samples = n1;             // select shortest file length
	else        nb_samples = n2;
	// nb_samples = 20000;
	printf("\n\nNumber of samples to be processed : %d \n",nb_samples);
	rewind(file_input_1);                    // reset the indices at the beginning of files
	rewind(file_input_2);

	// read input signal files
	input_1 = (double *)calloc(nb_samples,sizeof(double));
  	input_2 = (double *)calloc(nb_samples,sizeof(double));
	fread(input_1,sizeof(double),nb_samples,file_input_1);
	fread(input_2,sizeof(double),nb_samples,file_input_2);

	// Reserve and reset memory for estimated signals
	output_1 = (double *)calloc(nb_samples,sizeof(double));
  	output_2 = (double *)calloc(nb_samples,sizeof(double));
  	// Perform sources separation
	//separation(input_1, input_2, output_1, output_2, nb_samples, time_delay);
	separation_v0(input_1, input_2, output_1, output_2, nb_samples, time_delay);
	// Save outputs in files
	fwrite(output_1,sizeof(double),nb_samples,file_output_1);
	fwrite(output_2,sizeof(double),nb_samples,file_output_2);

	// Open and read reference signal files
	file_reference_1 = fopen("sobi_ref_1.txt","rb");
  	file_reference_2 = fopen("sobi_ref_2.txt","rb");
  	reference_1 = (double *)calloc(nb_samples,sizeof(double));
  	reference_2 = (double *)calloc(nb_samples,sizeof(double));
	fread(reference_1,sizeof(double),nb_samples,file_reference_1);
	fread(reference_2,sizeof(double),nb_samples,file_reference_2);

	// Evaluate estimation error
	cos_1   =  dot_product(reference_1,output_1,nb_samples)/norm(reference_1,nb_samples)/norm(output_1,nb_samples);
  	cos_2   =  dot_product(reference_2,output_2,nb_samples)/norm(reference_2,nb_samples)/norm(output_2,nb_samples);

  	//cos_1   =  dot_product(reference_1,output_2,nb_samples)/norm(reference_1,nb_samples)/norm(output_2,nb_samples);
  	//cos_2   =  dot_product(reference_2,output_1,nb_samples)/norm(reference_2,nb_samples)/norm(output_1,nb_samples);
	err_1   = 10*(log(1-(cos_1*cos_1))/log(10));
  	err_2   = 10*(log(1-(cos_2*cos_2))/log(10));
	printf("\nError index 1 (dB) = %f\n",     err_1);
	printf(  "Error index 2 (dB) = %f\n",     err_2);
	printf(  "cos(input_1, output_1) = %f\n",   cos_1);
	printf(  "cos(input_2, output_2) = %f\n\n", cos_2);

  // Free memory
  free(input_1);
  free(input_2);
  free(output_1);
  free(output_2);
  free(reference_1);
  free(reference_2);

  // Close files
  fclose(file_input_1);
  fclose(file_input_2);
  fclose(file_output_1);
  fclose(file_output_2);
  fclose(file_reference_1);
  fclose(file_reference_2);

  //printf("****************************************************************************\n");
  //printf("* For graphic results, see matlab program AffichageResultats.m sous Matlab *\n");
  //printf("****************************************************************************\n");

  return EXIT_SUCCESS;
}
