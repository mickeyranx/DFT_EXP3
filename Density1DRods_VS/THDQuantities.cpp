#include "THDQuantities.h"
#include <cmath>
#include <iostream>

double THDQuantities::psizeroD(double x)
{
	if (x >= 1) {
		x = 0.9999999;
	}
	return x + (1 - x) * log(1 - x);
}



double THDQuantities::psiZeroDPrime(double x)
{
	if (x >= 1) {
		x = 0.9999999;
	}
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
	//provides values from V to M - V - L - 1 
	for (int i = V; i < M - V - L + 3 ; i++)
	{
		//assign current rho
		double rho = rho_s[i];
		if(rho == 0){
			sum += psizeroD(weightedDensityOne(L, i, rho_s)) - psizeroD(weightedDensityZero(L, i, rho_s));

		}
		else {
			double a = rho * (log(rho) - 1)
				+ psizeroD(weightedDensityOne(L, i, rho_s))
				- psizeroD(weightedDensityZero(L, i, rho_s))
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
		sum += rho_s[i] - rho_0;
	}
	return sum;
}



double THDQuantities::weightedDensityOne(int L, int S, std::vector<double> &density_profile)
{
	double rho_s = 0;
	//sum over densities of L preceding lattice points
	for (int i = S - (L - 1); i < S + 1; i++)
	{
		rho_s += density_profile[i];
	}
	return rho_s;
}

double THDQuantities::weightedDensityZero(int L ,int S, std::vector<double> &density_profile)
{
	double rho_s = 0;
	//sum over densities of L-1 preceding lattice points
	for (int i = S - (L - 1); i < S; i++)
	{
		rho_s += density_profile[i];
	}

	return rho_s;
}



//maybe put muExS in different file
double THDQuantities::muExS(int L, int S, std::vector<double> &density_profile)
{
	double sum1 = 0;
	for (int i = S; i < S + L; i++)
	{
		sum1 += psiZeroDPrime(weightedDensityOne(L, i, density_profile));
	}

	double sum2 = 0;
	for (int i = S + 1; i < S + L; i++)
	{
		sum2 += psiZeroDPrime(weightedDensityZero(L, i, density_profile));
	}

	return sum1 - sum2;
}



