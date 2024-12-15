#include <stdio.h>
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

#include "../include/lbm_params.hpp"
#include "../include/function.hpp"
#include "../include/Ibm_3D_params.hpp"

int main(){

    // Velocity distribution function
    std::vector<double> f(Nx * Ny * 9, 0.0);

    // Density
    std::vector<double> rho(Nx * Ny, 1.0);

    // Velocity in x-direction
    std::vector<double> ux(Nx * Ny, 0.0);

    // Velocity in y-direction
    std::vector<double> uy(Nx * Ny, 0.0);

    // Local equilibrium distribution function
    std::vector<double> f_eq(Nx * Ny * 9, 0.0);

    double cu, etime;
    double uSqr;
    struct timeval t1, t2;

    for (int i = 0; i < Nx; i++){
        ux[idx_D2Q9(i, Ny - 1, Nx)] = U;
        for (int j = 0; j < Ny; j++){
            for (int k = 0; k < 9; k++){
                cu = c_D2Q9[k][0] * ux[idx_D2Q9(i, j, Nx)] + c_D2Q9[k][1] * uy[idx_D2Q9(i, j, Nx)];
                uSqr = ux[idx_D2Q9(i, j, Nx)] * ux[idx_D2Q9(i, j, Nx)] + uy[idx_D2Q9(i, j, Nx)] * uy[idx_D2Q9(i, j, Nx)];
                f[idx_D2Q9(i, j, k, Nx)] = w_D2Q9[k] * rho[idx_D2Q9(i, j, Nx)] * (1.0 + cu / c_s_square + cu * cu / (2 * c_s_square * c_s_square) - uSqr / (2 * c_s_square));
            }
        }
    }

    // gettimeofday(&t1, NULL);
    // for (int i = 0; i < 10000; i++){
    //     D2Q9_serial(f, rho, ux, uy, f_eq);
    // }
    // gettimeofday(&t2, NULL);
    // etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
    // etime = etime / 1000;
    // printf("Serial done, took %f sec.", etime);
   



    gettimeofday(&t1, NULL);
    for (int i = 0; i < 10000; i++){
        D2Q9_parallel(f, rho, ux, uy, f_eq);
    }
    gettimeofday(&t2, NULL);
    etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
    etime = etime / 1000;
    printf("Parallel done, took %f sec.\n", etime);

    /*
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            for (int k = 0; k < 9; k++){
                std::cout << f[idx_D2Q9(i, j, k, Nx)] << " ";
            }
            std::cout << "/ ";
        }
        std::cout << std::endl;
    }*/

    saveToCSV_D2Q9("lbm_results.csv", Nx, Ny, rho, ux, uy);



    // // Velocity distribution function
    // std::vector<double> f_3D(D3Q19::Nx * D3Q19::Ny * D3Q19::Nz * D3Q19::Q, 0.0);
    // // Density
    // std::vector<double> rho_3D(D3Q19::Nx * D3Q19::Ny * D3Q19::Nz, 1.0);
    // // Velocity in x-direction
    // std::vector<double> ux_3D(D3Q19::Nx * D3Q19::Ny * D3Q19::Nz, 0.0);
    // // Velocity in y-direction
    // std::vector<double> uy_3D(D3Q19::Nx * D3Q19::Ny * D3Q19::Nz, 0.0);
    // // Velocity in z-direction
    // std::vector<double> uz_3D(D3Q19::Nx * D3Q19::Ny * D3Q19::Nz, 0.0);
    // // Local equilibrium distribution function
    // std::vector<double> f_eq_3D(D3Q19::Nx * D3Q19::Ny * D3Q19::Nz * D3Q19::Q, 0.0);

    // initialize_D3Q19(f_3D, rho_3D, ux_3D, uy_3D, uz_3D);
    // gettimeofday(&t1, NULL);
    // for (int t = 0; t < 1000; t++) {
    //     D3Q19_parallel(f_3D, rho_3D, ux_3D, uy_3D, uz_3D, f_eq_3D);
    //     //std::cout << "Step " << t << " completed." << std::endl;
    // }
    // gettimeofday(&t2, NULL);
    // etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
    // etime = etime / 1000;
    // printf("Parallel_3D done, took %f sec.\n", etime);
    // saveToCSV_D3Q19("lbm_results_3D.csv", rho_3D, ux_3D, uy_3D, uz_3D);

    return 0;
}
