#include <fstream>
#include <filesystem>
#include <iostream>
#include <omp.h>
#include <string>

#include "../include/LBM_2D.hpp"

// Constructor
LBM_2D::LBM_2D(int Nx, int Ny, int Re, double U, double rho0)
    : Nx(Nx), Ny(Ny), Re(Re), U(U), rho0(rho0),
      f(Nx * Ny * 9, 0.0), f_temp(Nx * Ny * 9, 0.0),
      rho(Nx * Ny, rho0), ux(Nx * Ny, 0.0), uy(Nx * Ny, 0.0) {}

void LBM_2D::initialize() {
    double cu, uSqr;
    for (int i = 0; i < Nx; i++) {
        ux[idx_D2Q9(i, Ny - 1, Nx)] = U;
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < 9; k++) {
                cu = c_D2Q9[k][0] * ux[idx_D2Q9(i, j, Nx)] + c_D2Q9[k][1] * uy[idx_D2Q9(i, j, Nx)];
                uSqr = ux[idx_D2Q9(i, j, Nx)] * ux[idx_D2Q9(i, j, Nx)] + uy[idx_D2Q9(i, j, Nx)] * uy[idx_D2Q9(i, j, Nx)];
                f[idx_D2Q9(i, j, k, Nx)] = w_D2Q9[k] * rho[idx_D2Q9(i, j, Nx)] * (1.0 + cu / c_s_square + cu * cu / (2 * c_s_square * c_s_square) - uSqr / (2 * c_s_square));
            }
        }
    }
}

void LBM_2D::D2Q9_serial_iterate(int step) {
    double cu;
    double uSqr;
    double f_eq;
    for (int s = 0; s < step; s++) {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                // Collision
                for (int k = 0; k < 9; k++) {
                    cu = c_D2Q9[k][0] * ux[idx_D2Q9(i, j, Nx)] + c_D2Q9[k][1] * uy[idx_D2Q9(i, j, Nx)];
                    uSqr = ux[idx_D2Q9(i, j, Nx)] * ux[idx_D2Q9(i, j, Nx)] + uy[idx_D2Q9(i, j, Nx)] * uy[idx_D2Q9(i, j, Nx)];
                    f_eq = w_D2Q9[k] * rho[idx_D2Q9(i, j, Nx)] * (1.0 + cu / c_s_square + cu * cu / (2 * c_s_square * c_s_square) - uSqr / (2 * c_s_square));
                    f_temp[idx_D2Q9(i, j, k, Nx)] = (1.0 - 1.0 / tau) * f[idx_D2Q9(i, j, k, Nx)] + (1.0 / tau) * f_eq;
                }
            }
        }

        // Streaming
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                for (int k = 0; k < 9; k++) {
                    if (i - c_D2Q9[k][0] >= 0 && j - c_D2Q9[k][1] >= 0 && i - c_D2Q9[k][0] < Nx && j - c_D2Q9[k][1] < Ny) {
                        f[idx_D2Q9(i, j, k, Nx)] = f_temp[idx_D2Q9(i - c_D2Q9[k][0], j - c_D2Q9[k][1], k, Nx)];
                    }
                }
            }
        }

        for (int j = 0; j < Ny - 1; j++) {
            // Left boundary
            f[idx_D2Q9(0, j, 5, Nx)] = f[idx_D2Q9(0, j, 7, Nx)];
            f[idx_D2Q9(0, j, 1, Nx)] = f[idx_D2Q9(0, j, 3, Nx)];
            f[idx_D2Q9(0, j, 8, Nx)] = f[idx_D2Q9(0, j, 6, Nx)];

            // Right boundary
            f[idx_D2Q9(Nx - 1, j, 6, Nx)] = f[idx_D2Q9(Nx - 1, j, 8, Nx)];
            f[idx_D2Q9(Nx - 1, j, 3, Nx)] = f[idx_D2Q9(Nx - 1, j, 1, Nx)];
            f[idx_D2Q9(Nx - 1, j, 7, Nx)] = f[idx_D2Q9(Nx - 1, j, 5, Nx)];
        }

        for (int i = 0; i < Nx; i++) {
            // Bottom boundary
            f[idx_D2Q9(i, 0, 5, Nx)] = f[idx_D2Q9(i, 0, 7, Nx)];
            f[idx_D2Q9(i, 0, 2, Nx)] = f[idx_D2Q9(i, 0, 4, Nx)];
            f[idx_D2Q9(i, 0, 6, Nx)] = f[idx_D2Q9(i, 0, 8, Nx)];
        }

        double density;
        for (int i = 0; i < Nx; i++) {
            density = f[idx_D2Q9(i, Ny - 1, 0, Nx)] + f[idx_D2Q9(i, Ny - 1, 1, Nx)] + f[idx_D2Q9(i, Ny - 1, 3, Nx)] + 2 * (f[idx_D2Q9(i, Ny - 1, 2, Nx)] + f[idx_D2Q9(i, Ny - 1, 5, Nx)] + f[idx_D2Q9(i, Ny - 1, 6, Nx)]);
            f[idx_D2Q9(i, Ny - 1, 4, Nx)] = f[idx_D2Q9(i, Ny - 1, 2, Nx)];
            f[idx_D2Q9(i, Ny - 1, 8, Nx)] = f[idx_D2Q9(i, Ny - 1, 6, Nx)] - 0.5 * (f[idx_D2Q9(i, Ny - 1, 1, Nx)] - f[idx_D2Q9(i, Ny - 1, 3, Nx)]) + 0.5 * density * U;
            f[idx_D2Q9(i, Ny - 1, 7, Nx)] = f[idx_D2Q9(i, Ny - 1, 5, Nx)] + 0.5 * (f[idx_D2Q9(i, Ny - 1, 1, Nx)] - f[idx_D2Q9(i, Ny - 1, 3, Nx)]) - 0.5 * density * U;
        }

        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                rho[idx_D2Q9(i, j, Nx)] = 0;
                ux[idx_D2Q9(i, j, Nx)] = 0;
                uy[idx_D2Q9(i, j, Nx)] = 0;
                for (int k = 0; k < 9; k++) {
                    rho[idx_D2Q9(i, j, Nx)] += f[idx_D2Q9(i, j, k, Nx)];
                    ux[idx_D2Q9(i, j, Nx)] += c_D2Q9[k][0] * f[idx_D2Q9(i, j, k, Nx)];
                    uy[idx_D2Q9(i, j, Nx)] += c_D2Q9[k][1] * f[idx_D2Q9(i, j, k, Nx)];
                }
                ux[idx_D2Q9(i, j, Nx)] /= rho[idx_D2Q9(i, j, Nx)];
                uy[idx_D2Q9(i, j, Nx)] /= rho[idx_D2Q9(i, j, Nx)];
            }
        }
    }
}

void LBM_2D::D2Q9_parallel_iterate(int step) {
    double cu;
    double uSqr;
    double f_eq;
    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);
    for (int s = 0; s < step; s++) {
        // Collision
        #pragma omp parallel for private(cu, uSqr, f_eq) schedule(static) collapse(2)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                for (int k = 0; k < 9; k++) {
                    cu = c_D2Q9[k][0] * ux[idx_D2Q9(i, j, Nx)] + c_D2Q9[k][1] * uy[idx_D2Q9(i, j, Nx)];
                    uSqr = ux[idx_D2Q9(i, j, Nx)] * ux[idx_D2Q9(i, j, Nx)] + uy[idx_D2Q9(i, j, Nx)] * uy[idx_D2Q9(i, j, Nx)];
                    f_eq = w_D2Q9[k] * rho[idx_D2Q9(i, j, Nx)] * (1.0 + cu / c_s_square + cu * cu / (2 * c_s_square * c_s_square) - uSqr / (2 * c_s_square));
                    f_temp[idx_D2Q9(i, j, k, Nx)] = (1.0 - 1.0 / tau) * f[idx_D2Q9(i, j, k, Nx)] + (1.0 / tau) * f_eq;
                }
            }
        }

        // Streaming
        #pragma omp parallel for schedule(static) collapse(2)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                for (int k = 0; k < 9; k++) {
                    if (i - c_D2Q9[k][0] >= 0 && j - c_D2Q9[k][1] >= 0 && i - c_D2Q9[k][0] < Nx && j - c_D2Q9[k][1] < Ny) {
                        f[idx_D2Q9(i, j, k, Nx)] = f_temp[idx_D2Q9(i - c_D2Q9[k][0], j - c_D2Q9[k][1], k, Nx)];
                    }
                }
            }
        }

        for (int j = 0; j < Ny - 1; j++) {

            // Left boundary
            f[idx_D2Q9(0, j, 5, Nx)] = f[idx_D2Q9(0, j, 7, Nx)];
            f[idx_D2Q9(0, j, 1, Nx)] = f[idx_D2Q9(0, j, 3, Nx)];
            f[idx_D2Q9(0, j, 8, Nx)] = f[idx_D2Q9(0, j, 6, Nx)];

            // Right boundary
            f[idx_D2Q9(Nx - 1, j, 6, Nx)] = f[idx_D2Q9(Nx - 1, j, 8, Nx)];
            f[idx_D2Q9(Nx - 1, j, 3, Nx)] = f[idx_D2Q9(Nx - 1, j, 1, Nx)];
            f[idx_D2Q9(Nx - 1, j, 7, Nx)] = f[idx_D2Q9(Nx - 1, j, 5, Nx)];
        }

        for (int i = 0; i < Nx; i++) {
            // Bottom boundary
            f[idx_D2Q9(i, 0, 5, Nx)] = f[idx_D2Q9(i, 0, 7, Nx)];
            f[idx_D2Q9(i, 0, 2, Nx)] = f[idx_D2Q9(i, 0, 4, Nx)];
            f[idx_D2Q9(i, 0, 6, Nx)] = f[idx_D2Q9(i, 0, 8, Nx)];
        }

        double density;
        #pragma omp parallel for private(density)
        for (int i = 0; i < Nx; i++) {
            density = f[idx_D2Q9(i, Ny - 1, 0, Nx)] + f[idx_D2Q9(i, Ny - 1, 1, Nx)] + f[idx_D2Q9(i, Ny - 1, 3, Nx)] + 2 * (f[idx_D2Q9(i, Ny - 1, 2, Nx)] + f[idx_D2Q9(i, Ny - 1, 5, Nx)] + f[idx_D2Q9(i, Ny - 1, 6, Nx)]);
            f[idx_D2Q9(i, Ny - 1, 4, Nx)] = f[idx_D2Q9(i, Ny - 1, 2, Nx)];
            f[idx_D2Q9(i, Ny - 1, 8, Nx)] = f[idx_D2Q9(i, Ny - 1, 6, Nx)] - 0.5 * (f[idx_D2Q9(i, Ny - 1, 1, Nx)] - f[idx_D2Q9(i, Ny - 1, 3, Nx)]) + 0.5 * density * U;
            f[idx_D2Q9(i, Ny - 1, 7, Nx)] = f[idx_D2Q9(i, Ny - 1, 5, Nx)] + 0.5 * (f[idx_D2Q9(i, Ny - 1, 1, Nx)] - f[idx_D2Q9(i, Ny - 1, 3, Nx)]) - 0.5 * density * U;
        }

        #pragma omp parallel for schedule(static) collapse(2)
        // Combine two nested loops to improve load balancing
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                rho[idx_D2Q9(i, j, Nx)] = 0;
                ux[idx_D2Q9(i, j, Nx)] = 0;
                uy[idx_D2Q9(i, j, Nx)] = 0;
                for (int k = 0; k < 9; k++) {
                    rho[idx_D2Q9(i, j, Nx)] += f[idx_D2Q9(i, j, k, Nx)];
                    ux[idx_D2Q9(i, j, Nx)] += c_D2Q9[k][0] * f[idx_D2Q9(i, j, k, Nx)];
                    uy[idx_D2Q9(i, j, Nx)] += c_D2Q9[k][1] * f[idx_D2Q9(i, j, k, Nx)];
                }
                ux[idx_D2Q9(i, j, Nx)] /= rho[idx_D2Q9(i, j, Nx)];
                uy[idx_D2Q9(i, j, Nx)] /= rho[idx_D2Q9(i, j, Nx)];
            }
        }
    }
}

void LBM_2D::save_to_CSV(const std::string &filename) const {
    std::filesystem::create_directories("result");
    std::ofstream file("result/" + filename);
    file << "x,y,rho,ux,uy\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = j * Nx + i;
            file << i << "," << j << "," << rho[idx] << "," << ux[idx] << "," << uy[idx] << "\n";
        }
    }
    file.close();
}