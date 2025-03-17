#pragma once
#include <vector>


class THDQuantities
{

public:

	 static double psizeroD(double x);

	 static double psiZeroDPrime(double x);

	 /* calculates the canonical potential of a homogene density profile
	 note that wall_pos_1 and 2 are actually the positions after and before the barriers respectively*/
	 static double canonialPotHomogene(double rho_0, int wall_pos_1, int wall_pos_2, int rod_length);

	static double canonicalPot(std::vector<double> rho_s, std::vector<double> ext_pot, double mu, int rod_length, int wall_pos, int lattice_length);
	
	 //calculates the canonical potential of a density profile with a given external potential and a chemical potential 
	
	 //calculates the pressure of a homogene density distribution
	 static double surfacePressure(double rho_0, int rod_length);

	 //calculates the exakt Adsorbtion for a 1D lattice with Rods with determined density distribution rho_s a starting density rho_0 and the 2 positions of the barriers
	 static double exactAdsorbtionOneDim(std::vector<double> rho_s, double rho_0, int wall_pos_1, int wall_pos_2);

	 static double weightedDensityOne(int L, int S, std::vector<double>& density_profile);

	 static double weightedDensityZero(int L, int S, std::vector<double>& density_profile);

	 static double muExS(int L, int S, std::vector<double>& density_profile);

	 


};





