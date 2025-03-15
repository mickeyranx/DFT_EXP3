#include "THDQuantities.h"
#include "Lattice.h"
#include <cmath>
#include <iostream>

double THDQuantities::psizeroD(double x)
{
	return x + (1 - x) * log(1 - x);
}



double THDQuantities::psiZeroDPrime(double x)
{
	return -1*log(1-x);
}


double THDQuantities::canonialPotHomogene(double rho_0, int wall_pos_1, int wall_pos_2, int rod_length)
{
	return -1 * surfacePressure(rho_0, rod_length) * (wall_pos_2 - wall_pos_1 + 1);
}

double THDQuantities::canonicalPot(std::vector<double> rho_s, std::vector<double> ext_pot, double mu, int rod_length,int wall_pos, int lattice_length){
	int M = lattice_length;
	int V = wall_pos;
	int L = rod_length;
	double sum = 0;
	//liefert erst ab V und bis M - V - L - 1 Beitr�ge
	for (int i = V; i < M - V - L + 3 ; i++)
	{
		//assign current rho
		double rho = rho_s[i];
		if(rho == 0){
			sum += psizeroD(Lattice::weightedDensityOne(i, rho_s)) - psizeroD(Lattice::weightedDensityZero(i, rho_s));

		}
		else {
			double a = rho * (log(rho) - 1)
				+ psizeroD(Lattice::weightedDensityOne(i, rho_s))
				- psizeroD(Lattice::weightedDensityZero(i, rho_s))
				- (mu) * rho;
			
			sum += a;
		}

	}


	return sum;

}

double THDQuantities::surfacePressure(double rho_0 , int rod_length) {
	int L = rod_length;
	return log((1 - ((L - 1) * rho_0)) / (1 - L * rho_0));
}

double THDQuantities::exactAdsorbtionOneDim(std::vector<double> rho_s, double rho_0,  int wall_pos_1, int wall_pos_2) {
	double sum = 0;
	for (int i = wall_pos_1 + 1; i < wall_pos_2; i++)
	{
		//std::cout << "i = " << i << std::endl;
		//std::cout << rho_s[i] << std::endl;
		sum += rho_s[i] - rho_0;
	}
	return sum;
}




