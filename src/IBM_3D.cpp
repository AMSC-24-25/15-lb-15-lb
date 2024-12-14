#include <vector>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/time.h>
#include "../include/Ibm_3D_params.hpp"


void initialize_D3Q19(std::vector<double>& f, std::vector<double>& rho,
                      std::vector<double>& ux, std::vector<double>& uy,
                      std::vector<double>& uz) {
    double cu, uSqr;

    // Initialising density and speed
    #pragma omp parallel for collapse(3)
    for (int z = 0; z < Nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {

                // Setting the speed of the top boundary
                if (z == Nz - 1) {
                    ux[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] = U;
                }

                for (int k = 0; k < Q; k++) {
                    cu = c_D3Q19[k][0] * ux[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] +
                         c_D3Q19[k][1] * uy[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] +
                         c_D3Q19[k][2] * uz[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)];
                    uSqr = ux[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] * ux[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] +
                           uy[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] * uy[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] +
                           uz[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] * uz[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)];
                    f[idx_D3Q19(x, y, z, k, Nx, Ny, Nz)] = w_D3Q19[k] * rho[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] *
                        (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uSqr);
                }
            }
        }
    }
}



// 碰撞与流动函数
void D3Q19_parallel(std::vector<double>& f, std::vector<double>& rho,
                    std::vector<double>& ux, std::vector<double>& uy,
                    std::vector<double>& uz, std::vector<double>& f_eq) {
    double cu, uSqr;
    std::vector<double> f_temp(Nx * Ny * Nz * Q, 0.0);

    // 碰撞阶段
    #pragma omp parallel for private(cu, uSqr) collapse(3)
    for (int z = 0; z < Nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                for (int k = 0; k < Q; k++) {
                    cu = c_D3Q19[k][0] * ux[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] +
                         c_D3Q19[k][1] * uy[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] +
                         c_D3Q19[k][2] * uz[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)];
                    uSqr = ux[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] * ux[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] +
                           uy[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] * uy[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] +
                           uz[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] * uz[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)];

                    f_eq[idx_D3Q19(x, y, z, k, Nx, Ny, Nz)] = w_D3Q19[k] * rho[idx_D3Q19(x, y, z, 0, Nx, Ny, Nz)] *
                        (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uSqr);
                    
                    f_temp[idx_D3Q19(x, y, z, k, Nx, Ny, Nz)] = (1.0 - 1.0 / tau) * 
                        f[idx_D3Q19(x, y, z, k, Nx, Ny, Nz)] + (1.0 / tau) * f_eq[idx_D3Q19(x, y, z, k, Nx, Ny, Nz)];
                }
            }
        }
    }

    // 流动阶段，镜像反射
    #pragma omp parallel for collapse(3)
    for (int z = 0; z < Nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                for (int k = 0; k < Q; k++) {
                    int x_next = x - c_D3Q19[k][0];
                    int y_next = y - c_D3Q19[k][1];
                    int z_next = z - c_D3Q19[k][2];
                    if (x_next >= 0 && x_next < Nx && y_next >= 0 && y_next < Ny && z_next >= 0 && z_next < Nz) {
                        f[idx_D3Q19(x, y, z, k, Nx, Ny, Nz)] = 
                            f_temp[idx_D3Q19(x_next, y_next, z_next, k, Nx, Ny, Nz)];
                    }
                }
            }
        }
    }

    // 顶边和底边
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
            int idx_top = idx_D3Q19(x, y, Nz - 1, 0, Nx, Ny, Nz);
            int idx_bottom = idx_D3Q19(x, y, 0, 0, Nx, Ny, Nz);

            // 顶边设置速度
            ux[idx_top] = U;
            uy[idx_top] = 0.0;
            uz[idx_top] = 0.0;

            // 底边直接反弹
            for (int k = 0; k < Q; k++) {
                if (c_D3Q19[k][2] < 0) { // 处理向下的速度分布
                    f[idx_bottom + k] = f[idx_bottom + (Q - k - 1)];
                }
            }
        }
    }

}

void saveToCSV(const std::string& filename, const std::vector<double>& rho,
               const std::vector<double>& ux, const std::vector<double>& uy,
               const std::vector<double>& uz) {
    std::ofstream file(filename);
    file << "x,y,z,rho,ux,uy,uz\n";
    for (int z = 0; z < Nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int idx = idx_D3Q19(x, y, z, 0, Nx, Ny, Nz);
                file << x << "," << y << "," << z << ","
                     << rho[idx] << "," << ux[idx] << ","
                     << uy[idx] << "," << uz[idx] << "\n";
            }
        }
    }
    file.close();
}
int main(){

    struct timeval t1, t2;
    double etime;

    
    // Velocity distribution function
    std::vector<double> f(Nx * Ny * Nz * Q, 0.0);
    // Density
    std::vector<double> rho(Nx * Ny * Nz, 1.0);
    // Velocity in x-direction
    std::vector<double> ux(Nx * Ny * Nz, 0.0);
    // Velocity in y-direction
    std::vector<double> uy(Nx * Ny * Nz, 0.0);
    // Velocity in z-direction
    std::vector<double> uz(Nx * Ny * Nz, 0.0);
    // Local equilibrium distribution function
    std::vector<double> f_eq(Nx * Ny * Nz * Q, 0.0);

    initialize_D3Q19(f, rho, ux, uy, uz);
    gettimeofday(&t1, NULL);
    for (int t = 0; t < 1000; t++) {
        D3Q19_parallel(f, rho, ux, uy, uz, f_eq);
        std::cout << "Step " << t << " completed." << std::endl;
    }
    gettimeofday(&t2, NULL);
    etime = (t2.tv_sec - t1.tv_sec) * 1000 + (t2.tv_usec - t1.tv_usec) / 1000;
    etime = etime / 1000;
    printf("Parallel_3D done, took %f sec.", etime);
    saveToCSV("lbm_results_3D.csv", rho, ux, uy, uz);
   



}