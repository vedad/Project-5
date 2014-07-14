## Simulating a cold uniform spherical collapse
Here I simulate the evolution of systems such as star clusters. The objective is
to solve the N-body problem; a system bound by gravitational
attraction. The objects start with zero velocity (cold start) and are
distributed uniformly in a sphere of radius R. The equations of motion are
solved numerically by applying the ODE solver Runge-Kutta 4.


Here follows a list describing the function of each program.
* `gaussiandeviate` draws random numbers that were used in determining the
	initial positions and the masses of the bodies.
* `makeObjects` uses `gaussiandeviate` to draw random positions inside a
	sphere of radius R. It draws masses with a given mean value and deviation.
* `CelestialObject` contains the attributes of each object.
* `SolarSystem` solves the equation of motion for the entire system by using
	Runge-Kutta 4. The name of the program sticked from an earlier
	project where we solved the equation of motion for the solar system.
* `readData, readKinPotData, readObjectsData` are Python scripts for plotting
	the data (positions, energies, etc.).
