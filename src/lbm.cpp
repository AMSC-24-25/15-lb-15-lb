#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include <vector>
#include <iostream>

#include "../include/lbm_params.hpp"

/*  Used to calculate the index of D2Q9 index
    eg: f[idx_D2Q9(x, y, i, Nx)] = 1.0; 
    (x, y) means the coordinate of grid points
    i means the velocity direction at (x, y)
    Nx means the number of grid points in the x-direction
*/
inline size_t idx_D2Q9(int x, int y, int i, int Nx){
    return static_cast<size_t>((y * Nx + x) * 9 + i);
}

/*  eg: rho[idx_D2Q9(x, y, Nx)] = 1.0;
*/
inline size_t idx_D2Q9(int x, int y, int Nx){
    return static_cast<size_t>((y * Nx + x) * 9);
}

int main(){

    // Velocity distribution function
    std::vector<double> f(Nx * Ny * 9, 0.0);

    // Density
    std::vector<double> rho(Nx * Ny, 0.0);

    // Velocity in x-direction
    std::vector<double> ux(Nx * Ny, 0.0);

    // Velocity in x-direction
    std::vector<double> uy(Nx * Ny, 0.0);

    // Local equilibrium distribution function
    std::vector<double> f_eq(Nx * Ny * 9, 0.0);

    return 0;
}