#include <iostream>
#include "SolarSystem.hpp"
#include <ctime>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {
	
	int N = 200;
	double dt = 0.001;
	double R0 = 20;
	double tMax = 5;
	clock_t start, finish;	
	start = clock();

	/*	For the two-body problem. */
//	SolarSystem mySystem = SolarSystem("../data/parameters/Sirius.dat");
//	mySystem.setCenterOfMassPosition();
//	mySystem.setTotalMomentum();
//	mySystem.systemSimulation(dt, tMax, "leapfrog");
	
	/*  For the cluster simulation. */
	SolarSystem myCluster = SolarSystem(N, R0);
	myCluster.systemSimulation(dt, tMax, "leapfrog");

	finish = clock();
	cout << "Computation time: " << double(finish - start)/CLOCKS_PER_SEC << " seconds" << endl;
	
	/* Writing parameters to an information file. */
	fstream outFile;
	outFile.open("../data/info.dat", ios::out);
	outFile << N << " " << R0 << " " << dt << " " << tMax << endl;
	outFile.close();
	cout << "Completed!" << endl;

	return 0;
}


