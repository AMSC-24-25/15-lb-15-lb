#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/time.h>

#include "../include/LBM_2D.hpp"
#include "../include/LBM_3D.hpp"

int main(){

    int Nx = 100;
    int Ny = 100;
    int Nz = 100;
    int Re = 100;
    double U = 0.5;
    double rho0 = 1.0;

    double etime;
    struct timeval t1, t2;

    LBM_2D lbm_parallel(Nx, Ny, Re, U, rho0);
    lbm_parallel.initialize();

    gettimeofday(&t1, NULL);
    lbm_parallel.D2Q9_parallel_iterate(10000);
    gettimeofday(&t2, NULL);
    etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
    etime = etime / 1000;
    std::cout << "D2Q9_parallel done, took " << etime << " sec." << std::endl;

    lbm_parallel.save_to_CSV("lbm_results.csv");

    U = 0.3;
    LBM_3D lbm_3D(Nx, Ny, Nz, U, rho0);
    lbm_3D.initialize();

    gettimeofday(&t1, NULL);
    lbm_3D.D3Q19_parallel_iterate(100);
    gettimeofday(&t2, NULL);
    etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
    etime = etime / 1000;
    std::cout << "D3Q19_parallel done, took " << etime << " sec." << std::endl;

    lbm_3D.save_to_CSV("lbm_results_3D.csv");

    return 0;
}
