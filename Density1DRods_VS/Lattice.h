#pragma once

#include <vector>
//works in c++ 11+
#include <tuple>

class Lattice
{
public:
	//rod length 
	static int L;
	//Lattice Length
	static int M; 
	//Wall position (2nd wall is at M - V)
	static int V;
	

	static double rho_s_0;

	static double mu_ex;


	//picard mixing ------------------------------
	//alpha for picard as a global var
	static double alpha;
	//maximum alpha for picard
	static const double ALPHA_MAX;
	//minimum alpha for picard
	static const double ALPHA_MIN;
	
	static double epsilon_old; 
	//-------------------------------------------
	static std::vector<double> wall_profile;
	
	

	//creates the 1D lattice (elements are densities) with a given size
	static std::vector<double> createLattice(int M);
	//creates the wall-profile of the lattice (elements are zero (outside walls) or 1 (within the walls))
	static std::vector<int> createWallProfile(int lattice_size, int wall_start);

	static double weightedDensityOne(int S, std::vector<double> density_profile);

	static double weightedDensityZero(int S, std::vector<double> density_profile);

	static double muExS(int S, std::vector<double> density_profile);

	//calculates next alpha with convergence condition
	static void setNextAlpha(double epsilon_new);

	//this is the actual main function of the project
	static std::vector<double> calculateDensityProfile(int rod_length, double nu_0, int lattice_size, int wall_position, int steps, double alpha);
	//calculates the next iteration
	static std::vector<double> iterate(std::vector<double> current_density_profile);

	static void checkAdsorbtionCondition(std::vector<double> rho_s, double delta_mu);

	//static std::tuple<std::vector<double>, double epsilon_old> 

	//return the density profile of the lattice after N iterations
	
	
	


};

