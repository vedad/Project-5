#include <iostream>
#include "SolarSystem.hpp"
#include <ctime>
#include <fstream>

using namespace std;

int main(int argc, char* argv[]) {
	
	clock_t start, finish;	
	start = clock();
//	SolarSystem mySystem = SolarSystem("../data/parameters/sunJupiter.dat");
	SolarSystem myCluster = SolarSystem(100, 20);
//	mySystem.setCenterOfMassPosition();
//	mySystem.setTotalMomentum();
	myCluster.systemSimulation(0.01, 5, true, false, "leapfrog");
	finish = clock();
	cout << "Computation time: " << double(finish - start)/CLOCKS_PER_SEC << " seconds" << endl;

	fstream outFile;
	outFile.open("../data/info.dat", ios::out);
	outFile << 100 << " " << 20 << " " << 0.01 << " " << 5 << endl;
	outFile.close();

	return 0;
}


