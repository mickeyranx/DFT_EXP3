// Density1DRods_VS.cpp : Diese Datei enthält die Funktion "main". Hier beginnt und endet die Ausführung des Programms.
//

#include <iostream>
#include "Lattice.h"
#include <fstream>
#include <iomanip>

int main()
{
	//rodlength
	int L = 3;
	//number of lattice points
	int M = 208;
	//wall position (2nd wall (M-V))
	int V = 10;
	//initial density
	double nu_0 = 0.6;
	//parameter for picard mixing
	double starting_alpha = 0.01;

	int manual_max_steps = 200;

	std::vector<double> rho = Lattice::calculateDensityProfile(L, nu_0, M, V, manual_max_steps, starting_alpha);
	
	//create textfile
	std::ofstream Data("data3.txt");
	//put data into textfile
	for (auto density : rho)
	{
		Data << std::setprecision(8) << density << "\n";
	}

	
	


	//TODO: daten exportieren
   
}


