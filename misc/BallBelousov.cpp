// BallBelousov.cpp : Ball (1994) model in 2D phase plane
// uses library built from rk4.cpp

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "rk4.h"

//Reaction rates
double k1 = 1.0;
double k2 = 1.0;
double k3 = 1.0;
//Reactant sum
double S = 1.0;

//Derivative functions
//x and y are functions of the three reactant concentrations and the sum of concentration
double xprime(std::vector<double> x, double time) {
    double xp = (k3 * (S - x[1] - x[0]) * (S + 2 * x[0])) - (k2 * (S + x[1] - x[0]) * (S + 2 * x[0]));
    return xp;
}

double yprime(std::vector<double> x, double time) {
    double yp = k2 * (S + x[1] - x[0]) * (S + 2 * x[0]) - 2 * k1 * (S - x[1] - x[0]) * (S + x[1] - x[0]) + k3 * (S - x[1] - x[0]) * (S + 2 * x[0]);
    return yp;
}

int main()
{
    int niter;
    int npoints = 1; // number of starting points
    double dx = 0.1; // size of step away from equilibrium point
    double dt;
    //Initial state
    double x0, y0;
    std::cout << "Enter k1" << std::endl;
    std::cin >> k1;
    std::cout << "Enter k2" << std::endl;
    std::cin >> k2;
    std::cout << "Enter k3" << std::endl;
    std::cin >> k3;
  
    std::cout << "Enter number of points" << std::endl;
    std::cin >> npoints;
    std::cout << "Enter point spread distance" << std::endl;
    std::cin >> dx;
    std::cout << "Enter number of steps" << std::endl;
    std::cin >> niter;
    std::cout << "Enter step size" << std::endl;
    std::cin >> dt;
    
    // equilibrium position
    x0 = S * (2 * k1 - k2 - k3) / (2 * (k1 + k2 + k3));
    y0 = 3 * S * (k3 - k2) / (2 * (k1 + k2 + k3));

    //Create vector of differential equations
    std::vector<double (*)(std::vector<double>, double)> diffeqs;
    diffeqs.push_back(xprime);
    diffeqs.push_back(yprime);

    // Solve from initial conditions
    for (int k = 1; k < npoints + 1; k++) {
        std::string fname = std::to_string(k);
        //Initialize solver
        std::vector<double> init{ x0 + k*dx, y0 + k*dx };
        std::cout << "Starting at " << init[0] << " " << init[1] << std::endl;
        rk4::rk4Solve Solver(init, 0.0, dt, diffeqs);

        //Output files
        std::fstream file;
        file.open(fname.append(".csv"), std::fstream::out);
        if (!file.is_open()) { std::cerr << "Failed to create output file." << std::endl; }

        for (int j = 0; j < niter; j++) {
            Solver.iterate();
            for (double x : Solver.stateCurrent) {
                std::cout << x << " ";
                file << x << ",";
            }
            file << "\n";
            std::cout << std::endl;
        }
        file.close();
    }
    
    return 0;
}
