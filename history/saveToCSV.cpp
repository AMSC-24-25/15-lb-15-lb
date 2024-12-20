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
    for (int z = 0; z < Nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int idx = (z * Ny + y) * Nx + x;
                file << x << "," << y << "," << z << ","
                     << rho[idx] << "," << ux[idx] << ","
                     << uy[idx] << "," << uz[idx] << "\n";
            }
        }
    }
    file.close();
}
