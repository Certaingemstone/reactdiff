// Runge-Kutta, 4th order
#pragma once
#include <vector>

namespace rk4
{
	class rk4Solve
	{
	public:
		//State information
		std::vector<double> stateCurrent;
		//Time
		double t, dt;
		//Number of state variables
		size_t n;
		//Pointers to functions taking in current state and time and outputting derivative for each state variable.
		std::vector<double (*)(std::vector<double>, double)> derivatives;
		//constructor
		rk4Solve( std::vector<double> state, double time, double timestep, std::vector<double (*)(std::vector<double>, double)> inFunctions);
		//Takes current state, calculates k1, k2, k3, k4, returns next state and updates object state variables.
		std::vector<double> iterate();
	};
}