/* ************************************************************************** */
/*                                                                            */
/*    Copyright 2013 Stefano Zamuner                                          */
/*                                                                            */
/*    This file is part of BioCpp.                                            */
/*                                                                            */
/*    BioCpp is free software: you can redistribute it and/or modify          */
/*    it under the terms of the GNU General Public License as published by    */
/*    the Free Software Foundation, either version 3 of the License, or       */
/*    (at your option) any later version.                                     */
/*                                                                            */
/*    BioCpp is distributed in the hope that it will be useful,               */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*    GNU General Public License for more details.                            */
/*                                                                            */
/*    You should have received a copy of the GNU General Public License       */
/*    along with BioCpp.  If not, see <http://www.gnu.org/licenses/>.         */
/*                                                                            */
/* ************************************************************************** */

#ifndef BIOCPP_RNG_H
#define BIOCPP_RNG_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string>

struct rngs{
  gsl_rng * rng;
  const gsl_rng_type * type;

  std::string seed_filename;

  rngs(std::string filename = "seed_file.dat", unsigned long int seed = 0 );

  int RandomUniformInteger(int min_val, int max_val);  /* values in [ min_val, max_val ) */
  double RandomUniformDouble(double min_val, double max_val);
  double RandomGaussianDouble(double average, double sigma);
  void SaveStatus();
};

inline rngs::rngs(std::string filename, unsigned long int seed){
  seed_filename = filename;
	type = gsl_rng_ranlxs2;
	rng = gsl_rng_alloc (type);
  gsl_rng_set ( rng, seed);
	FILE * seed_file_in;
	seed_file_in = fopen(seed_filename.c_str(), "r");
	if (seed_file_in == NULL){
		FILE *seed_file_out;
		seed_file_out = fopen(seed_filename.c_str(), "w");
		gsl_rng_fwrite (seed_file_out, rng);
		fclose(seed_file_out);
	}
	seed_file_in = fopen(seed_filename.c_str(), "r");
	gsl_rng_fread (seed_file_in, rng);
	fclose(seed_file_in);
}

inline int rngs::RandomUniformInteger( int min_val,int max_val ){ 
	double ran = gsl_rng_uniform( rng )*(  max_val - min_val + 1);
	return  (int)ran + min_val;
}

inline double rngs::RandomUniformDouble(double min_val, double max_val){
	double ran = gsl_rng_uniform( rng )*(  max_val - min_val );
	return ran + min_val;
}

inline double rngs::RandomGaussianDouble(double average, double sigma){
  double ran = gsl_ran_gaussian(rng, sigma);
	return ran + average;
}

inline void rngs::SaveStatus(){
	FILE * seed_file_out;
	seed_file_out = fopen(seed_filename.c_str(), "w");
	gsl_rng_fwrite (seed_file_out, rng);
	fclose(seed_file_out);
}

#endif
