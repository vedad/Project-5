#include <iostream>
#include <armadillo>
#include <cmath>
#include "makeObjects.hpp"
#include "gaussiandeviate.hpp"

using namespace std;
using namespace arma;

const double PI = 3.1415926535;

vec createPosition(long seed, double R) {

	double u = ran2(&seed);
	double v = ran2(&seed);
	double w = ran2(&seed);

	double phi = 2.*PI*w;
	double theta = acos(1-2*v);
	double r = R * pow(u,0.333);
	
	vec a = zeros<vec>(3);
	a[0] = r*sin(theta)*cos(phi);
	a[1] = r*sin(theta)*sin(phi);
	a[2] = r*cos(theta);

	return a;
}

vec createVelocity() {

	vec b = zeros<vec>(3);
	return b;
}

double createMass(long seed, double mean, double deviation) {

	return mean + deviation * gaussian_deviate(&seed);
}
