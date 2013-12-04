#include <iostream>
#include <armadillo>
#include "CelestialObject.hpp"

using namespace std;
using namespace arma;

const double EPSILON = 0.15;

CelestialObject :: CelestialObject(string name, vec position, vec velocity, double mass) {

	this->mass = mass;
	this->position = position;
	this->velocity = velocity;
	this->name = name;

}
/*
CelestialObject :: CelestialObject(vec position, vec velocity, double mass) {

	this->mass = mass;
	this->position = position;
	this->velocity = velocity;
}
*/

// Set-functions
void CelestialObject :: setVelocity(vec newVelocity) {

	velocity = newVelocity;

}

void CelestialObject :: setPosition(vec newPosition) {
	position = newPosition;
}


// Get-functions
vec CelestialObject :: getDistanceTo(CelestialObject other) {
	
	vec distance = other.getPosition() - this->getPosition();
	return distance;

}

vec CelestialObject :: getForce(CelestialObject other) {
	
	vec r = this->getDistanceTo(other);
	double absR = norm(r,2);
	vec force = (this->getMass() * other.getMass() / (pow(absR,2) + pow(EPSILON,2))) * r / absR;
	return force;
		
}

vec CelestialObject :: getAcceleration(CelestialObject other) {

	vec force = getForce(other);
	vec acceleration = force/(this->getMass());
	return acceleration;

}
   
double CelestialObject :: getKineticEnergy(CelestialObject object) {
	double kineticEnergy = 0.5 * object.getMass() * norm(object.getVelocity(),2) * norm(object.getVelocity(),2);
	return kineticEnergy;
}

 
double CelestialObject :: getPotentialEnergy(CelestialObject other) {

	vec r = this->getDistanceTo(other);
	double potentialEnergy = - this->getMass() * other.getMass() / norm(r,2); 
	return potentialEnergy;
}

vec CelestialObject :: getVelocity() { return velocity; }
vec CelestialObject :: getPosition() { return position; }
double CelestialObject :: getMass() { return mass; }
string CelestialObject :: getName() { return name; }

