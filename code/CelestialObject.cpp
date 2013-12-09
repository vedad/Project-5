#include <iostream>
#include <armadillo>
#include "CelestialObject.hpp"

using namespace std;
using namespace arma;

const double EPSILON_SQ = 0.0225;

CelestialObject :: CelestialObject(string name, vec position, vec velocity, double mass) {

	this->mass = mass;
	this->position = position;
	this->velocity = velocity;
	this->name = name;

}

/* Set-functions */

// Sets a new velocity for an object.
void CelestialObject :: setVelocity(vec newVelocity) {

	velocity = newVelocity;

}

// Sets a new position for an object.
void CelestialObject :: setPosition(vec newPosition) {
	position = newPosition;
}


/* Get-functions */

// Gets the distance between two objects.
vec CelestialObject :: getDistanceTo(CelestialObject other) {
	
	vec distance = other.getPosition() - this->getPosition();
	return distance;

}

// Gets the force between two objects.
vec CelestialObject :: getForce(CelestialObject other) {
	
	vec r = this->getDistanceTo(other);
	double absR = norm(r,2);
	vec force = (this->getMass() * other.getMass() / (pow(absR,2) + EPSILON_SQ)) * r / absR;
	return force;
		
}

// Gets the acceleration of an object given the force from another.
vec CelestialObject :: getAcceleration(CelestialObject other) {

	vec force = getForce(other);
	vec acceleration = force/(this->getMass());
	return acceleration;

}

// Gets the kinetic energy of an object.
double CelestialObject :: getKineticEnergy(CelestialObject object) {
	double kineticEnergy = 0.5 * object.getMass() * norm(object.getVelocity(),2) * norm(object.getVelocity(),2);
	return kineticEnergy;
}

// Gets the potential energy between two objects. 
double CelestialObject :: getPotentialEnergy(CelestialObject other) {

	vec r = this->getDistanceTo(other);
	double potentialEnergy = - this->getMass() * other.getMass() / norm(r,2); 
	return potentialEnergy;
}

vec CelestialObject :: getVelocity() { return velocity; }
vec CelestialObject :: getPosition() { return position; }
double CelestialObject :: getMass() { return mass; }
string CelestialObject :: getName() { return name; }

