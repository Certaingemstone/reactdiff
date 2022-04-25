// Runge-Kutta, 4th order

#include "rk4.h"

namespace rk4
{
	rk4Solve::rk4Solve(std::vector<double> state, double time, double timestep, std::vector<double (*)(std::vector<double>, double)> inFunctions) 
		: t(time), dt(timestep), n(state.size()), stateCurrent(state), derivatives(inFunctions) { }

	std::vector<double> rk4Solve::iterate()
	{
		std::vector<double> stateNext (n), stateIntermediate (n);
		std::vector<double> k1 (n), k2 (n), k3 (n), k4 (n);

		//Obtain k1 (derivative vector at current state)
		for (int i = 0; i < n; i++) {
			k1[i] = derivatives[i](stateCurrent, t);
		}

		//Obtain k2 (derivative vector at half timestep, linear from stateCurrent along k1) 
		for (int i = 0; i < n; i++) {
			stateIntermediate[i] = stateCurrent[i] + k1[i] * dt / 2.0;
		}
		for (int i = 0; i < n; i++) {
			k2[i] = derivatives[i](stateIntermediate, t + dt/2.0);
		}

		//Obtain k3 (derivative vector at half timestep, linear from stateCurrent along k2)
		for (int i = 0; i < n; i++) {
			stateIntermediate[i] = stateCurrent[i] + k2[i] * dt / 2.0;
		}
		for (int i = 0; i < n; i++) {
			k3[i] = derivatives[i](stateIntermediate, t + dt/2.0);
		}

		//Obtain k4 (derivative vector at full timestep, linear from stateCurrent along k3)
		for (int i = 0; i < n; i++) {
			stateIntermediate[i] = stateCurrent[i] + k3[i] * dt;
		}
		for (int i = 0; i < n; i++) {
			k4[i] = derivatives[i](stateIntermediate, t + dt);
		}

		//Update stateNext using k
		for (int i = 0; i < n; i++) {
			stateNext[i] = stateCurrent[i] + ((k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0) * dt;
		}
		//Update stateCurrent and time
		stateCurrent = stateNext;
		t = t + dt;
		return stateNext;
	}
}