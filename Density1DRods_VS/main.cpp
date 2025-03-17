// Density1DRods_VS.cpp : Diese Datei enthält die Funktion "main". Hier beginnt und endet die Ausführung des Programms.
//

#include "THDQuantities.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <tuple>



//-----------------------------------
//          global params
//-----------------------------------
double const STARTING_ALPHA = 0.01;
double const ALPHA_MIN = 0.0001;
double const ALPHA_MAX = 0.03;
//Wall position 
int const V = 10;
//call count of startSimulation() to manage data writing
static int callcount = 0;
static int callcount2 = 0;

//return the next alpha for the next iteration
static double nextAlpha(double epsilon_old, double epsilon_new, double alpha) {
	//do nothing if this is the first iteration (alpha stays)
	if (epsilon_old == -1) {	
		return alpha;
	}
	else {
		//if new epsilon is smaller than old multiply by 1.1 otherwise divide by 5
		(epsilon_new < epsilon_old) ? (alpha = 1.1 * alpha) : (alpha = alpha / 5);
		
		if (alpha > ALPHA_MAX) {
			return 0.03;
		}
		else if (alpha < ALPHA_MIN) {
			return 0.0001;
		
		}
		else {
			return alpha;
		}
	}
}

static std::tuple<double, double> checkAdsorbtionCondition(std::vector<double> rho_s, double delta_rho, double rho_s_0, int M, int L)
{
	//calculate and print the exat Adsorbtion
	// V - 1 is the actual position of the first wall and ((M - V) - L + 1) - 1 of the second
	//std::cout << "exakt Adsorption G = " << THDQuantities::exactAdsorbtionOneDim(rho_s, rho_s_0, V - 1, M - V - L) / 2 << std::endl;
	double G_ex = THDQuantities::exactAdsorbtionOneDim(rho_s, rho_s_0, V - 1, M - V - L) / 2;
	//estimate the Adsorption with a variation of rho_0 and integrating with respect to mu
	double mu_id_p = log(rho_s_0 + delta_rho);
	double mu_ex_p = -L * log(1 - (rho_s_0 + delta_rho) * L) + (L - 1) * log(1 - (rho_s_0 + delta_rho) * (L - 1));
	double mu_id_m = log(rho_s_0 + delta_rho);
	double mu_ex_m = -L * log(1 - (rho_s_0 - delta_rho) * L) + (L - 1) * log(1 - (rho_s_0 - delta_rho) * (L - 1));
	double mu_p = mu_id_p + mu_ex_p;
	double mu_m = mu_id_m + mu_ex_m;


	double gamma_p = (mu_ex_p - (2 * L - 1) * THDQuantities::surfacePressure(rho_s_0 + delta_rho, L)) / 2;
	double gamma_m = (mu_ex_m - (2 * L - 1) * THDQuantities::surfacePressure(rho_s_0 - delta_rho, L)) / 2;
	//double Gamma = THDQuantitie;
	//std::cout << "estimated Adsorbtion G = " << -1 * (gamma_p - gamma_m) / (mu_p - mu_m) << std::endl;
	double G_est = -1 * (gamma_p - gamma_m) / (mu_p - mu_m);
	return std::make_tuple(G_ex, G_est);
}




//returns the density profile and epsilon_new (needed to determine next alpha) of the next iteration
static std::tuple<std::vector<double>, double> iterate(
	std::vector<double>& wall,
	std::vector<double>& rho_prev,
	double alpha, double rho_s_0, double mu_ex, double epsilon_old ,int M, int L) {
	//---------------------------------
	//  calculate next density profile
	//---------------------------------
	std::vector<double> rho_new(M, 0);
	for (int i = L; i < M - L; i++)
	{
		rho_new[i] = rho_s_0 * exp(mu_ex - THDQuantities::muExS(L, i, rho_prev)) * wall[i];
	}

	//---------------------------------
	//         picard mixing
	//---------------------------------
	std::vector<double> rho_next(M, 0);
	double epsilon_new = 0;
	for (int i = 0; i < M; i++)
	
	{
		//prefix c = current
		double c_rho_new = rho_new[i];

		double c_rho = rho_prev[i];

		rho_next[i] = (1 - alpha) * c_rho + (alpha)*c_rho_new;

		epsilon_new += pow((c_rho_new - c_rho), 2) * 1 / M;


	}
	
	return std::make_tuple(rho_next, epsilon_new);

}





//this function initializes important lists, the first density profile and values and starts to iterate
static void startSimulation(int M ,int L , double nu_0, int max_iter, std::string filename1, std::string filename2) {

	//------------------------------------
	//           setup
	//------------------------------------
	//create wall profile (0 outside walls, 1 inside walls)
	std::vector<double> wall(M, 0);
	for (int i = 0; i < M; i++)
	{
		if ((i > V - 1) && (i < (M - V) - L)) { //if lattice point within walls 
			wall[i] = 1;
		}
		else {
			wall[i] = 0;
		}
	}

	//calculate starting values for starting density profile
	double rho_s_0 = nu_0 / L; 
	double mu_ex = -L * log(1 - rho_s_0 * L) + (L - 1) * log(1 - rho_s_0 * (L - 1));

	//calculate first (0th) density_profile
	std::vector<double> rho_0(M);
	for (int i = 0; i < M; i++)
	{
		rho_0[i] = rho_s_0 * wall[i];
	}

	//initialize values for first iteration
	double alpha = STARTING_ALPHA;
	double epsilon_old = -1; //set first epsilon to -1 to detect first iteration
	double epsilon_new = 0;
	std::vector<double> rho_prev = rho_0;
	std::vector<double> rho_next = rho_0;
	//------------------------------------
	//           iterations
	//------------------------------------
	for (int i = 0; i < max_iter; i++)
	{
		if (epsilon_old <= pow(10, -9) && epsilon_old > 0) break; //stop iterating if epsilon is very small

		//calculate next density profile and epsilon new
		std::tie(rho_next, epsilon_new) = iterate(wall, rho_prev, alpha, rho_s_0, mu_ex, epsilon_old, M, L); //unpack tuple
		//prepare for the next iteration
		rho_prev = rho_next;
		alpha = nextAlpha(epsilon_old, epsilon_new, alpha);
		epsilon_old = epsilon_new;
		
	}
	
	//------------------------------------
	//        write data to file
	//------------------------------------
	
	//------------------------------------
	//        1. density profile
	//------------------------------------


	if (callcount == 0) { //first function call
		std::ofstream file(filename1); //create new file
		file << std::fixed << std::setprecision(8);
		for (double density : rho_next) //put data into textfile
		{
			file << density << "\n";
		}
	}else{

		std::ifstream inFile(filename1);
		std::vector<std::string> lines(M, "");

		std::string line; 
		int row = 0;

		while (std::getline(inFile, line) && row < M) {
			lines[row] = line; 
			row++;
		}
		inFile.close();
		for (int i = 0; i < M; i++) //put data into textfile
		{
			std::stringstream ss; 
			ss << std::fixed << std::setprecision(8) << (rho_next[i]);
			if (!lines[i].empty()) {
				lines[i] += "\t";
			}
			lines[i] += ss.str();
		}

		std::ofstream outFile(filename1);
		for (const auto& l : lines)
		{
			outFile << l << std::endl;
		}
		outFile.close();
	}
	callcount++;

	//------------------------------------
	//    2. other quantities
	//------------------------------------
	if (callcount2 == 0) {
		std::ofstream file2(filename2);
		file2.close();
	}
	
	std::ofstream file2(filename2, std::ios::app); //open textfile in append mode
	file2 << "----------------------------------" << "\n";
	file2 << "     nu_0 = " << nu_0 << "          " << "\n";
	file2 << "----------------------------------" << "\n";
	file2 << std::fixed << std::setprecision(5);
	double delta_mu = 0.001;
	double mu = log(rho_s_0) + mu_ex;
	file2 << "mu =             " << mu << "\n";
	double omega_rho = THDQuantities::canonicalPot(rho_prev, wall, mu, L, V, M);
	double omega_rho_0 = THDQuantities::canonialPotHomogene(rho_s_0, V, M - V - L - 1, L);
	double gamma = THDQuantities::canonicalPot(rho_prev, wall, mu, L, V - 1, M) - THDQuantities::canonialPotHomogene(rho_s_0, V, M - V - L - 1, L);
	//double gamma =  THDQuantities::canonialPotHomogene(rho_s_0, V - 1, M - V - L, L);
	file2 << "pressure:        " << THDQuantities::surfacePressure(rho_s_0, L) << "\n";
	file2 << "omega_rho =      " << omega_rho << "\n";
	file2 << "omega_rho_0 =    " << omega_rho_0 << "\n";
	file2 << "gamma_ana =      " << (mu_ex - (2 * L - 1) * THDQuantities::surfacePressure(rho_s_0, L)) / 2 << "\n";
	file2 << "gamma =          " << (omega_rho - omega_rho_0) / 2 << "\n";
	
	double G_ex;
	double G_est;
	std::tie(G_ex, G_est) = checkAdsorbtionCondition(rho_prev, delta_mu, rho_s_0, M, L);
	file2 << "exakt Gamma:     " << G_ex << "\n";
	file2 << "estimated Gamma: " << G_est << "\n";
	callcount2++;

}






int main()
{
	clock_t start = clock();
	//--------------------------------
	//  set parameters for simulation
	//--------------------------------
	//rodlength
	int L = 3;
	//number of lattice points
	int M = 208;
	//wall position (2nd wall (M-V))
	int V = 10;
	//parameter for picard mixing
	double starting_alpha = 0.01;
	//max no. iterations
	int manual_max_steps = 5000;

	//--------------------------------
	//  simulate with given parameters
	//    from nu_0 = 0.1 to 0.9
	//--------------------------------
	//initial density
	double nu_0 = 0.1;
	for (int i = 0; i < 9; i++)
	{
		
		startSimulation(M, L, nu_0, manual_max_steps, "L" + std::to_string(L) + "_rho.txt", "L" + std::to_string(L) + "_gamma.txt");
		nu_0 += 0.1;
	}
	
	clock_t end = clock();
	double elapsed = double(end - start) / CLOCKS_PER_SEC;
	printf("exectution time: %.3f sec", elapsed);



}


