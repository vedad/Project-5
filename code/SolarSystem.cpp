#include <iostream>
#include <armadillo>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include "makeObjects.hpp"
#include "SolarSystem.hpp"

using namespace std;
using namespace arma;

const double DIMENSION = 3;
const double PI = 3.1415926535;
/* Need to change gravitational constant so that it can be used with AU as the
 * unit of length, and years as the unit of time */

SolarSystem :: SolarSystem(string systemfile) {

	fstream inFile;
	inFile.open("../data/parameters/sunJupiter.dat", ios::in);
	
	double x0, y0, z0, v0x, v0y, v0z, m;
	string name;

	while (inFile >> name >> x0 >> y0 >> z0 >> v0x >> v0y >> v0z >> m) {
		
		vec position; position << x0 << y0 << z0;
		vec velocity; velocity << v0x << v0y << v0z;
		CelestialObject newObject = CelestialObject(name, position, velocity, m);
		addObject(newObject);
		cout << "Added: " << newObject.getName() << endl;
	}
}
SolarSystem :: SolarSystem(int N, double R) {
	
	this->R = R;
	this->N = N;

	long seed1 = -1;
	long seed2 = -2;

	for (int i=0; i < N; i++) {
		
		vec velocity = createVelocity();
		vec position = createPosition(&seed1, R);
		double mass = createMass(&seed2, 10.0, 1.0);
		CelestialObject newObject = CelestialObject(position, velocity, mass);
		addObject(newObject);
	}

}

void SolarSystem :: addObject(CelestialObject newObject) {
	objects.push_back(newObject);
}

/* Advancing the system a time step dt by the Leap Frog algorithm */

void SolarSystem :: leapFrog(double dt) {
	
	mat newVelocity = zeros<mat>(DIMENSION, getNoOfObjects());
	mat newPosition = zeros<mat>(DIMENSION, getNoOfObjects());

	for (int i=0; i < getNoOfObjects(); i++) {
		newVelocity.col(i) = objects[i].getVelocity() + 0.5 * dt * getSystemAcceleration(objects[i]);
		newPosition.col(i) = objects[i].getPosition() + dt * newVelocity.col(i);
	}

	for (int i=0; i < getNoOfObjects(); i++) {
		objects[i].setVelocity(newVelocity.col(i));
		objects[i].setPosition(newPosition.col(i));
	}

	for (int i=0; i < getNoOfObjects(); i++) {
		newVelocity.col(i) = objects[i].getVelocity() + 0.5 * dt * getSystemAcceleration(objects[i]);
		objects[i].setVelocity(newVelocity.col(i));
	}
	
}

/* Advancing the system a time step dt by the RK4 method */

void SolarSystem :: advance(double dt) {

	mat velK1 = zeros<mat>(DIMENSION, getNoOfObjects());
	mat velK2 = zeros<mat>(DIMENSION, getNoOfObjects());
	mat velK3 = zeros<mat>(DIMENSION, getNoOfObjects());
	mat velK4 = zeros<mat>(DIMENSION, getNoOfObjects());

	mat posK1 = zeros<mat>(DIMENSION, getNoOfObjects());
	mat posK2 = zeros<mat>(DIMENSION, getNoOfObjects());
	mat posK3 = zeros<mat>(DIMENSION, getNoOfObjects());
	mat posK4 = zeros<mat>(DIMENSION, getNoOfObjects());

	mat thisVelocity = zeros<mat>(DIMENSION, getNoOfObjects());
	mat thisPosition = zeros<mat>(DIMENSION, getNoOfObjects());

	mat newVelocity = zeros<mat>(DIMENSION, getNoOfObjects());
	mat newPosition = zeros<mat>(DIMENSION, getNoOfObjects());

	double dt2 = dt / 2.0;
	double dt6 = dt / 6.0;
		
	for (int i=0; i < getNoOfObjects(); i++) {
		thisPosition.col(i) = objects[i].getPosition();
		thisVelocity.col(i) = objects[i].getVelocity();

		velK1.col(i) = getSystemAcceleration(objects[i]);
		posK1.col(i) = objects[i].getVelocity();			
	}
	
	// Move entire system by half the step length.
	for (int i=0; i < getNoOfObjects(); i++) {
		objects[i].setPosition(thisPosition.col(i) + dt2 * posK1.col(i));
	}
	
	// Calculate K2.
	for (int i=0; i < getNoOfObjects(); i++) {
		velK2.col(i) = getSystemAcceleration(objects[i]);
		posK2.col(i) = posK1.col(i) + dt2 * velK1.col(i);
	}
	
	// Move system by half a step length.
	for (int i=0; i < getNoOfObjects(); i++) {
		objects[i].setPosition(thisPosition.col(i) + dt2 * posK2.col(i));
	}
	
	// Calculate K3.
	for (int i=0; i < getNoOfObjects(); i++) {
		velK3.col(i) = getSystemAcceleration(objects[i]);
		posK3.col(i) = posK1.col(i) + dt2 * velK2.col(i);
	}
	
	// Move system by the step length.
	for (int i=0; i < getNoOfObjects(); i++) {
		objects[i].setPosition(thisPosition.col(i) + dt*posK3.col(i));
	}
	
	// Calculate K4.
	for (int i=0; i < getNoOfObjects(); i++) {
		velK4.col(i) = getSystemAcceleration(objects[i]);
		posK4.col(i) = posK1.col(i) + dt * velK3.col(i);
	}
	// Calculate the new velocity and position after a time step dt.
	for (int i=0; i < getNoOfObjects(); i++) {
		newVelocity.col(i) = thisVelocity.col(i) + dt6 * (velK1.col(i) + 2*velK2.col(i) + 2*velK3.col(i) + velK4.col(i));	
		newPosition.col(i) = thisPosition.col(i) + dt6 * (posK1.col(i) + 2*posK2.col(i) + 2*posK3.col(i) + posK4.col(i));

		objects[i].setVelocity(newVelocity.col(i));
		objects[i].setPosition(newPosition.col(i));
	}

}

void SolarSystem :: systemSimulation(double dt, double tMax, bool energy, bool angMom, string solver) {
	
	ofstream *newEnergyFile;
	vector<ofstream*> objectEnergyList;
	ofstream *newAngMomFile;
	vector<ofstream*> objectAngMomList;

	if (energy) {
		
		for (int i=0; i < getNoOfObjects(); i++) {
			
			ostringstream energyFile;
			energyFile << "../data/conservations/energy/" << objects[i].getName() << ".dat";
			newEnergyFile = new ofstream(energyFile.str().c_str());
//			*newEnergyFile << "Energy for: " << objects[i].getName() << endl;
			objectEnergyList.push_back(newEnergyFile);

		}
	}

	if (angMom) {

		for (int i=0; i < getNoOfObjects(); i++) {
			
			ostringstream angMomFile;
			angMomFile << "../data/conservations/angularmomentum/" << objects[i].getName() << ".dat";
			newAngMomFile = new ofstream(angMomFile.str().c_str());
			*newAngMomFile << "Angular momentum for: " << objects[i].getName() << endl;
			objectAngMomList.push_back(newAngMomFile);

		}
	}

	ofstream *newPositionFile;
	vector<ofstream*> objectFileList;

	for (int i=0; i < getNoOfObjects(); i++) {
		
		ostringstream ss;
		ss << i;
		string name = ss.str();
		ostringstream objectFile;
		objectFile << "../data/objects/" << name << ".dat";
		newPositionFile = new ofstream(objectFile.str().c_str());
		*newPositionFile << "Positions for: " << name << endl;
		objectFileList.push_back(newPositionFile);
				
	}

	cout << "Advancing system..." << endl;
	
	if (solver.compare("RK4") == 0) {

		double totalTime = 0;
		clock_t start, finish;

		for (double t=0; t < tMax; t+=dt) {
			
			start = clock();
			this->advance(dt);
			finish = clock();
			totalTime += double(finish - start)/CLOCKS_PER_SEC;
		
			for (int i=0; i < getNoOfObjects(); i++) {
				*objectFileList[i] << objects[i].getPosition()[0] << " " << objects[i].getPosition()[1] << " " << objects[i].getPosition()[2] << endl;

				if (energy) {
					*objectEnergyList[i] << t << " " << objects[i].getKineticEnergy(objects[i]) << " " << getSystemPotentialEnergy(objects[i]) << " " << getTotalEnergy(objects[i]) << endl;
				}
				if (angMom) {
					*objectAngMomList[i] << t << " " << getAngularMomentum(objects[i]) << endl;
				}
			}
		}
		double avgTimeStep = totalTime / (tMax / dt);
		cout << "Computation time for one timestep using RK4: " << avgTimeStep << " seconds" << endl;
	}

	if (solver.compare("leapfrog") == 0) {

		double totalTime = 0;
		clock_t start, finish;

		for (double t=0; t < tMax; t+=dt) {
		
			start = clock();
			this->leapFrog(dt);
			finish = clock();
			totalTime += double(finish - start)/CLOCKS_PER_SEC;
		
			for (int i=0; i < getNoOfObjects(); i++) {
				*objectFileList[i] << objects[i].getPosition()[0] << " " << objects[i].getPosition()[1] << " " << objects[i].getPosition()[2] << endl;

				if (energy) {
					*objectEnergyList[i] << t << " " << objects[i].getKineticEnergy(objects[i]) << " " << getSystemPotentialEnergy(objects[i]) << " " << getTotalEnergy(objects[i]) << endl;
				}
				if (angMom) {
					*objectAngMomList[i] << t << " " << getAngularMomentum(objects[i]) << endl;
				}
			}
		}
		double avgTimeStep = totalTime / (tMax / dt);
		cout << "Computation time for one timestep using Leapfrog: " << avgTimeStep << " seconds" << endl;
	}

	for (int i = 0; i < getNoOfObjects(); i++) {
		objectFileList[i]->close();
		if (energy) { objectEnergyList[i]->close(); }
		if (angMom) { objectAngMomList[i]->close(); }
	}
	
//	outFile.close();
}

// Calculating the force from all celestial objects on an object
vec SolarSystem :: getSystemForce(CelestialObject object) {
	
	// Initialize systemForce as a vector of size 2, with zeros as entries since
	// I use += function below.
	vec systemForce = zeros<vec>(DIMENSION);
	for (int i=0; i < getNoOfObjects(); i++) {
	
		if (object.getName() == objects[i].getName()) { continue; }

		else {

			vec r = object.getDistanceTo(objects[i]);
			systemForce += object.getForce(objects[i]);
		}
	}

	return (getGravConst() * systemForce);
}

// Calculating the acceleration of an object.
vec SolarSystem :: getSystemAcceleration(CelestialObject object) {
	
	vec systemForce = getSystemForce(object);
	vec systemAcceleration = systemForce / object.getMass();

	return systemAcceleration;
}

double SolarSystem :: getAngularMomentum(CelestialObject object) {

//	c v = object.getVelocity();
//	doubler = object.getDistanceTo(Sun);
//	double angularVelocity = norm(object.getVelocity(),2) / norm(object.getDistanceTo(Sun),2);
	vec v = object.getVelocity();
	vec r = object.getDistanceTo(objects[0]);
	double m = object.getMass();
	double detRV = r(0) * v(1) - r(1) * v(0);
	double angularMomentum = m * detRV;

	return angularMomentum;

}

// Calculating the total energy of an object.
double SolarSystem :: getTotalEnergy(CelestialObject object) {

	double totalEnergy = 0.0; //zeros<vec>(getNoOfObjects());

	double kineticEnergy = object.getKineticEnergy(object);
	double potentialEnergy = getSystemPotentialEnergy(object);

	totalEnergy = kineticEnergy + potentialEnergy;
	return totalEnergy;
}
/* 
vec SolarSystem :: getSystemKineticEnergy(CelestialObject object) {
	
	vec kineticEnergy = zeros<vec>(getNoOfObjects());
	for (int i=0; i < getNoOfObjects(); i++) {
		kineticEnergy(i) = getKineticEnergy(objects[i]);	
	}
	return kineticEnergy;
}
*/

double SolarSystem :: getSystemPotentialEnergy(CelestialObject object) {

	double potentialEnergy = 0.0; //zeros<vec>(getNoOfObjects());

	for (int i=0; i < getNoOfObjects(); i++) {
	
		if (object.getName() == objects[i].getName()) { continue; }

		else {
			vec r = object.getDistanceTo(objects[i]);
			potentialEnergy += object.getPotentialEnergy(objects[i]);
		}
	}
	return potentialEnergy;
}

vec SolarSystem :: getCenterOfMassPosition() {

	vec CM = zeros<vec>(DIMENSION);
	double totalMass = 0;

	for (int i=0; i < getNoOfObjects(); i++) {
		CM += objects[i].getMass() * objects[i].getPosition();
		totalMass += objects[i].getMass();
	}
	return CM/totalMass;
}

void SolarSystem :: setCenterOfMassPosition() {
	
	vec CM = getCenterOfMassPosition();

	for (int i=0; i < getNoOfObjects(); i++) {
		objects[i].setPosition(objects[i].getPosition() - CM);
	}
}

vec SolarSystem :: getTotalMomentum() {

	vec totalMomentum = zeros<vec>(DIMENSION);
	for (int i=0; i < getNoOfObjects(); i++) {
		totalMomentum += objects[i].getMass() * objects[i].getVelocity();
	}
	return totalMomentum;
}

void SolarSystem :: setTotalMomentum() {

//	objects[0].setVelocity(objects[0].getVelocity() * 0);
	vec totalMomentum = getTotalMomentum();
	objects[0].setVelocity(objects[0].getVelocity() - totalMomentum/objects[0].getMass());
}

double SolarSystem :: getGravConst() {
	
	totalMass = N * 10	
	double volume = (4./3) * PI * pow(R,3);
	double density = totalMass / volume;
	double gravConst = (3 * PI) / (32 * density);
	return gravConst;
}


int SolarSystem :: getNoOfObjects() { return objects.size(); }


