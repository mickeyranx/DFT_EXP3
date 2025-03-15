#include "Lattice.h"
#include <iostream>
#include "THDQuantities.h"
#include <string>
std::vector<double> Lattice::createLattice(int M)
{
    return std::vector<double>();
}

std::vector<int> Lattice::createWallProfile(int lattice_size, int wall_start)
{
    return std::vector<int>();
}

//default values (not used)
int Lattice::M = 0;
int Lattice::L = 0;
int Lattice::V = 0;

double Lattice::rho_s_0 = 0;
double Lattice::mu_ex = 0;

double Lattice::alpha = 0;

std::vector<double> Lattice::wall_profile = {};

//initialize alpha max and min 
const double Lattice::ALPHA_MAX = 0.03;

const double Lattice::ALPHA_MIN = 0.0001;
//set epsilon old to -1 to detect first iteration
double Lattice::epsilon_old = -1;

double Lattice::weightedDensityOne(int S, std::vector<double> density_profile)
{
    //TODO: check boundary cond?
    double rho_s = 0;
    //sum over densities of L preceding lattice points
    for (int i = S - (L-1); i < S+1; i++)
    {
        rho_s += density_profile[i];
    }
    return rho_s;
}

double Lattice::weightedDensityZero(int S, std::vector<double> density_profile)
{
    //TODO: check boundary cond?
    double rho_s = 0;
    //sum over densities of L-1 preceding lattice points
    for (int i = S - (L - 1); i < S; i++)
    {
        rho_s += density_profile[i];
    }
    
    return rho_s;
}

double Lattice::muExS(int S,  std::vector<double> density_profile)
{

    //sum
    // !!! check boundaries !!!
    double sum1 = 0;
    for (int i = S; i < S+L; i++)
    {
        sum1 += THDQuantities::psiZeroDPrime(weightedDensityOne(i, density_profile));
    }

    double sum2 = 0;
    for (int i = S+1; i < S+L; i++)
    {
        sum2 += THDQuantities::psiZeroDPrime(weightedDensityZero(i, density_profile));
    }

    return sum1 - sum2;
}

void Lattice::setNextAlpha(double epsilon_new) {
    //do nothing if this is the first iteration (alpha stays)
    if (epsilon_old == -1) {
        epsilon_old = epsilon_new;
        return;
    }
    else {
        //if new epsilon is smaller than old multiply by 1.1 otherwise divide by 5
        (epsilon_new < epsilon_old) ? (alpha = 1.1 * alpha) : (alpha = alpha / 5);
        epsilon_old = epsilon_new;
        if (alpha > ALPHA_MAX) {
            alpha = ALPHA_MAX;
            return;
        }
        else if (alpha < ALPHA_MIN) {
            alpha = ALPHA_MIN;
            return;
        }
        else {
            return;
        }
    }
}



std::vector<double> Lattice::calculateDensityProfile(int rod_length, double nu_0, int lattice_size, int wall_position , int steps, double starting_alpha)
{
    //assign parameters to global variables
    M = lattice_size;
    L = rod_length;
    V = wall_position;
    //note that V - 1 is the actual position of the first wall and ((M - V) - L + 1) - 1 of the second
    
    //create wallprofile and assign to global variable wallprofile
    std::vector<double> wall(M,0);
    for (int i = 0; i < M; i++)
    {

        //if lattice point within walls
        if((i > V - 1) && (i < (M-V)-L)){
            wall[i] = 1;
        }
        else {
            wall[i] = 0;
        }
    }
    wall_profile = wall;

    alpha = starting_alpha;
    

    //starting values
    rho_s_0 = nu_0 / L;
    mu_ex = -L * log(1 - rho_s_0 * L) + (L - 1) * log(1 - rho_s_0 * (L - 1));

    //calculate first (0th) density_profile
    std::vector<double> rho_0(M);
    for (int i = 0; i < M; i++)
    {
        rho_0[i] = rho_s_0 * wall_profile[i];
    }


    std::vector<double> rho = rho_0;
    std::vector<double> rho_prev = rho_0;
    for (int i = 0; i < steps; i++)
    {
        //std::cout << "iteration = " << i << ", epsilion = " << epsilon_old << std::endl;
        //stop iteration if epsilon is very small
        if (epsilon_old <= pow(10, -9) && epsilon_old > 0) break;
        rho = iterate(rho_prev);
        rho_prev = rho;
    }

    //check for adsorbtion 
    //TODO: 
    //calculate pressure 
    double delta_mu = 0.001;
    double mu = log(rho_s_0) + mu_ex;
    std::cout << "mu = " << mu << std::endl;
    double omega_rho = THDQuantities::canonicalPot(rho, wall_profile, mu, L, V, M);
    double omega_rho_0 = THDQuantities::canonialPotHomogene(rho_s_0, V, M - V - L - 1, L);
    double gamma = THDQuantities::canonicalPot(rho, wall_profile,mu , L, V - 1, M) - THDQuantities::canonialPotHomogene(rho_s_0, V, M - V - L - 1, L);
    //double gamma =  THDQuantities::canonialPotHomogene(rho_s_0, V - 1, M - V - L, L);
    std::cout << "pressure = " << THDQuantities::surfacePressure(rho_s_0, L) << std::endl;
    
    std::cout << "gamma_ana = " << (mu_ex - (2 * L - 1) * THDQuantities::surfacePressure(rho_s_0, L))/2 << std::endl;
    std::cout << "omega_rho = " << omega_rho << std::endl;
    std::cout << "omega_rho_0 = " << omega_rho_0 << std::endl;
    std::cout << "gamma = " << (omega_rho - omega_rho_0)/2 << std::endl;
    checkAdsorbtionCondition(rho, delta_mu);



    return rho;
}

std::vector<double> Lattice::iterate(std::vector<double> current_density_profile)
{
    //assign current density profile to a shorter variable
    std::vector<double> rho = current_density_profile;
    //calculate rho new
    std::vector<double> rho_new(M);
    // !!! Check for Boundaries !!!
    for (int i = L; i < M - L; i++)
    {
        rho_new[i] = rho_s_0 * exp(mu_ex - muExS(i, rho)) * wall_profile[i];
    }

    std::vector<double> rho_next(M);
    double epsilon_new = 0;
    //picard mixing
    for (int i = 0; i < M; i++)
    {
        //prefix c = current
        double c_rho_new = rho_new[i];
        
        double c_rho = rho[i];

        rho_next[i] = (1-alpha) * c_rho + (alpha) * c_rho_new;

        epsilon_new += pow((c_rho_new - c_rho), 2) * 1 / M;

    }

    //determine and set next alpha 
    setNextAlpha(epsilon_new);

    //return newly calculated profile
    return rho_next;

}

void Lattice::checkAdsorbtionCondition(std::vector<double> rho_s, double delta_rho)
{

    //calculate and print the exat Adsorbtion
    // V - 1 is the actual position of the first wall and ((M - V) - L + 1) - 1 of the second
    std::cout << "exakt Adsorption G = " << THDQuantities::exactAdsorbtionOneDim(rho_s, rho_s_0, V - 1, M - V - L) / 2 << std::endl;

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
    std::cout << "estimated Adsorbtion G = " <<  -1* (gamma_p - gamma_m)/(mu_p - mu_m) << std::endl;

}

