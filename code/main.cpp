#include <iostream>
#include "SolarSystem.hpp"
#include <ctime>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {
	
	int N = 100;
	double dt = 0.001;
	double R0 = 20;
	double tMax = 5;
	clock_t start, finish;	
	start = clock();
//	SolarSystem mySystem = SolarSystem("../data/parameters/sunJupiter.dat");
	SolarSystem myCluster = SolarSystem(N, R0);
//	mySystem.setCenterOfMassPosition();
//	mySystem.setTotalMomentum();
	myCluster.systemSimulation(dt, tMax, "leapfrog");
//	mySystem.systemSimulation(dt, tMax, "leapfrog");
	finish = clock();
	cout << "Computation time: " << double(finish - start)/CLOCKS_PER_SEC << " seconds" << endl;
	
	fstream outFile;
	outFile.open("../data/info.dat", ios::out);
	outFile << N << " " << R0 << " " << dt << " " << tMax << endl;
	outFile.close();
	cout << "Still here." << endl;

	return 0;
}


