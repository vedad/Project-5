#include <iostream>
#include "SolarSystem.hpp"
#include <ctime>

using namespace std;

int main(int argc, char* argv[]) {
	
	clock_t start, finish;	
	start = clock();
//	SolarSystem mySystem = SolarSystem("../data/parameters/sunJupiter.dat");
	SolarSystem myCluster = SolarSystem(100, 20);
//	mySystem.setCenterOfMassPosition();
//	mySystem.setTotalMomentum();
	myCluster.systemSimulation(0.01, 2, false, false, "RK4");
	finish = clock();
	cout << "Computation time: " << double(finish - start)/CLOCKS_PER_SEC << " seconds" << endl;

	return 0;
}


