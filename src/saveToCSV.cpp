#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include "../include/Ibm_3D_params.hpp"


void saveToCSV_D2Q9(const std::string& filename, int Nx, int Ny, 
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

void saveToCSV_D3Q19(const std::string& filename, const std::vector<double>& rho,
               const std::vector<double>& ux, const std::vector<double>& uy,
               const std::vector<double>& uz) {
    std::ofstream file(filename);
    file << "x,y,z,rho,ux,uy,uz\n";
    for (int z = 0; z < D3Q19::Nz; z++) {
        for (int y = 0; y < D3Q19::Ny; y++) {
            for (int x = 0; x < D3Q19::Nx; x++) {
                int idx = idx_D3Q19(x, y, z, 0, D3Q19::Nx, D3Q19::Ny, D3Q19::Nz);
                file << x << "," << y << "," << z << ","
                     << rho[idx] << "," << ux[idx] << ","
                     << uy[idx] << "," << uz[idx] << "\n";
            }
        }
    }
    file.close();
}