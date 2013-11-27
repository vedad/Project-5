#include <armadillo>

using namespace std;
using namespace arma;

class CelestialObject {

	public:
		CelestialObject(string, vec, vec, double);
		string name;
		vec velocity, position;
		double mass;

		// Set-functions
		void setVelocity(vec);
		void setPosition(vec);
		
		// Get-functions
		vec getDistanceTo(CelestialObject);
		vec getForce(CelestialObject);
		vec getAcceleration(CelestialObject);
		vec getVelocity();
		vec getPosition();
		string getName();
		double getMass();
		double getKineticEnergy(CelestialObject);
		double getPotentialEnergy(CelestialObject);

};

