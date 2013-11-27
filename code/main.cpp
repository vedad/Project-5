#include <iostream>
#include "SolarSystem.hpp"
#include <ctime>

using namespace std;

int main(int argc, char* argv[]) {
	
	clock_t start, finish;	
	start = clock();
	SolarSystem mySystem = SolarSystem("../data/parameters/sunJupiter.dat");
	mySystem.setCenterOfMassPosition();
	mySystem.setTotalMomentum();
	mySystem.systemSimulation(1e-2, 500, false, false, "leapfrog");
	finish = clock();
	cout << "Computation time: " << double(finish - start)/CLOCKS_PER_SEC << " seconds" << endl;

	return 0;
}


