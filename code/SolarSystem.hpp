#include <vector>
#include "CelestialObject.hpp"

using namespace std;
using namespace arma;

class SolarSystem {
	
	public:
		SolarSystem(string);
		SolarSystem(int, double);
		vec getSystemForce(CelestialObject);
		vec getSystemAcceleration(CelestialObject);
		void addObject(CelestialObject);
		void advance(double);
		void leapFrog(double);
		void systemSimulation(double, double, bool, bool, string);
		double getSystemPotentialEnergy(CelestialObject);
		double getTotalEnergy(CelestialObject);
		double getAngularMomentum(CelestialObject);
		int getNoOfObjects();
		vec getCenterOfMassPosition();
		void setCenterOfMassPosition();
		vec getTotalMomentum();
		void setTotalMomentum();
		double getGravConst();
		int getBoundObjects(CelestialObject);
		double getClusterEnergy();
		double getBoundClusterEnergy();
		double getClusterPotentialEnergy();
		double getClusterKineticEnergy();

	
		vector<CelestialObject> objects;
		double R;
		int N;

};
