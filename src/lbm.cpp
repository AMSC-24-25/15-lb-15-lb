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

void saveToCSV(const std::string& filename, int Nx, int Ny, 
               const std::vector<double>& rho, 
               const std::vector<double>& ux, 
               const std::vector<double>& uy) {
    std::ofstream file(filename);
    file << "x,y,rho,ux,uy\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = j * Nx + i;
            file << i << "," << j << "," << rho[idx] << "," << ux[idx] << "," << uy[idx] << "\n";
        }
    }
    file.close();
}

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
                f[idx_D2Q9(i, j, k, Nx)] = w_D2Q9[k] * rho[idx_D2Q9(i, j, Nx)] * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uSqr);
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
    printf("Parallel done, took %f sec.", etime);

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

    saveToCSV("lbm_results.csv", Nx, Ny, rho, ux, uy);

    return 0;
}
