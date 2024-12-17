

void initializeData(std::vector<double>& f, std::vector<double>& rho, std::vector<double>& ux, std::vector<double>& uy, std::vector<double>& uz) {
    double cu, uSqr;
    for (int z = 0; z < Nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int idx_rho = idx_D3Q19(x,y,z,Nx,Ny);
                int idx_f = idx_D3Q19(x,y,z,0, Nx,Ny);
                
                rho[idx_rho] = 1.0;  // 初始密度
                ux[idx_rho] = 0.0;   // 初始速度
                uy[idx_rho] = 0.0;
                uz[idx_rho] = 0.0;

                if (z == Nz - 1) {  // 顶部边界条件
                    ux[idx_rho] = U;
                }

                for (int k = 0; k < 19; k++) {
                    cu = c_D3Q19[k][0] * ux[idx_rho] + c_D3Q19[k][1] * uy[idx_rho] + c_D3Q19[k][2] * uz[idx_rho];
                    uSqr = ux[idx_rho] * ux[idx_rho] + uy[idx_rho] * uy[idx_rho] + uz[idx_rho] * uz[idx_rho];
                    f[idx_f + k] = w_D3Q19[k] * rho[idx_rho] * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uSqr);
                }
            }
        }
    }
}



void D3Q19_parallel(std::vector<double>& f, std::vector<double>& rho, std::vector<double>& ux, std::vector<double>& uy, std::vector<double>& uz, std::vector<double>& f_eq){

    double cu, uSqr;
    std::vector<double> f_temp(Nx * Ny * Nz * 19, 0.0);
    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);


    // Collision
    #pragma omp parallel for private(cu, uSqr) schedule(static) collapse(3)
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            for (int z = 0; z < Nz; z++){
                for (int k = 0; k < 19; k++){
                    cu = c_D3Q19[k][0] * ux[idx_D3Q19(i, j, z, Nx, Ny)] +
                         c_D3Q19[k][1] * uy[idx_D3Q19(i, j, z, Nx, Ny)] +
                         c_D3Q19[k][2] * uz[idx_D3Q19(i, j, z, Nx, Ny)];
                    uSqr = ux[idx_D3Q19(i, j, z, Nx, Ny)] * ux[idx_D3Q19(i, j, z, Nx, Ny)] +
                           uy[idx_D3Q19(i, j, z, Nx, Ny)] * uy[idx_D3Q19(i, j, z, Nx, Ny)] +
                           uz[idx_D3Q19(i, j, z, Nx, Ny)] * uz[idx_D3Q19(i, j, z, Nx, Ny)];
                    f_eq[idx_D3Q19(i, j, z, k, Nx, Ny)] = w_D3Q19[k] * rho[idx_D3Q19(i, j, z, Nx, Ny)] * 
                        (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uSqr);
                    f_temp[idx_D3Q19(i, j, z, k, Nx, Ny)] = (1.0 - 1.0 / tau) * f[idx_D3Q19(i, j, z, k, Nx, Ny)] + 
                                                       (1.0 / tau) * f_eq[idx_D3Q19(i, j, z, k, Nx, Ny)];
                }
            }
        }
    }

    // Streaming
    #pragma omp parallel for schedule(static) collapse(3)
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            for (int z = 0; z < Nz; z++){
                for (int k = 0; k < 19; k++){
                    int i_new = i - c_D3Q19[k][0];
                    int j_new = j - c_D3Q19[k][1];
                    int z_new = z - c_D3Q19[k][2];
                    if(i_new >= 0 && j_new >= 0 && z_new >= 0 && 
                       i_new < Nx && j_new < Ny && z_new < Nz){
                        f[idx_D3Q19(i, j, z, k, Nx, Ny)] = std::max(0.0,f_temp[idx_D3Q19(i_new, j_new, z_new, k, Nx, Ny)]);
                    }
                }
            } 
        }
    }

   

    // 左右边界条件
    for (int j = 0; j < Ny; j++) {
        for (int z = 0; z < Nz; z++) {
            // 左边界 (x=0)
            f[idx_D3Q19(0, j, z, 1, Nx, Ny)] = f[idx_D3Q19(0, j, z, 2, Nx, Ny)];
            f[idx_D3Q19(0, j, z, 5, Nx, Ny)] = f[idx_D3Q19(0, j, z, 6, Nx, Ny)];
            f[idx_D3Q19(0, j, z, 8, Nx, Ny)] = f[idx_D3Q19(0, j, z, 9, Nx, Ny)];

            // 右边界 (x=Nx-1)
            f[idx_D3Q19(Nx - 1, j, z, 2, Nx, Ny)] = f[idx_D3Q19(Nx - 1, j, z, 1, Nx, Ny)];
            f[idx_D3Q19(Nx - 1, j, z, 6, Nx, Ny)] = f[idx_D3Q19(Nx - 1, j, z, 5, Nx, Ny)];
            f[idx_D3Q19(Nx - 1, j, z, 9, Nx, Ny)] = f[idx_D3Q19(Nx - 1, j, z, 8, Nx, Ny)];
        }
    }

    // 前后边界条件
    for (int i = 0; i < Nx; i++) {
        for (int z = 0; z < Nz; z++) {
            // 前边界 (y=0)
            f[idx_D3Q19(i, 0, z, 3, Nx, Ny)] = f[idx_D3Q19(i, 0, z, 4, Nx, Ny)];
            f[idx_D3Q19(i, 0, z, 7, Nx, Ny)] = f[idx_D3Q19(i, 0, z, 10, Nx, Ny)];
            f[idx_D3Q19(i, 0, z, 11, Nx, Ny)] = f[idx_D3Q19(i, 0, z, 12, Nx, Ny)];

            // 后边界 (y=Ny-1)
            f[idx_D3Q19(i, Ny - 1, z, 4, Nx, Ny)] = f[idx_D3Q19(i, Ny - 1, z, 3, Nx, Ny)];
            f[idx_D3Q19(i, Ny - 1, z, 10, Nx, Ny)] = f[idx_D3Q19(i, Ny - 1, z, 7, Nx, Ny)];
            f[idx_D3Q19(i, Ny - 1, z, 12, Nx, Ny)] = f[idx_D3Q19(i, Ny - 1, z, 11, Nx, Ny)];
        }
    }

    // 底部边界条件
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            f[idx_D3Q19(i, j, 0, 5, Nx, Ny)] = f[idx_D3Q19(i, j, 0, 6, Nx, Ny)];
            f[idx_D3Q19(i, j, 0, 13, Nx, Ny)] = f[idx_D3Q19(i, j, 0, 14, Nx, Ny)];
            f[idx_D3Q19(i, j, 0, 15, Nx, Ny)] = f[idx_D3Q19(i, j, 0, 16, Nx, Ny)];
        }
    }
    double density;

    // 顶部边界条件 (速度 U)
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            // 计算密度
            density = f[idx_D3Q19(i, j, Nz - 1, 0, Nx, Ny)] +
                      f[idx_D3Q19(i, j, Nz - 1, 1, Nx, Ny)] + f[idx_D3Q19(i, j, Nz - 1, 2, Nx, Ny)] +
                      f[idx_D3Q19(i, j, Nz - 1, 3, Nx, Ny)] + f[idx_D3Q19(i, j, Nz - 1, 4, Nx, Ny)] +
                      2 * (f[idx_D3Q19(i, j, Nz - 1, 5, Nx, Ny)] + f[idx_D3Q19(i, j, Nz - 1, 6, Nx, Ny)]);

            // 设置顶部边界的分布函数
            f[idx_D3Q19(i, j, Nz - 1, 6, Nx, Ny)] = f[idx_D3Q19(i, j, Nz - 1, 5, Nx, Ny)];
            f[idx_D3Q19(i, j, Nz - 1, 14, Nx, Ny)] = f[idx_D3Q19(i, j, Nz - 1, 13, Nx, Ny)] - density * U;
            f[idx_D3Q19(i, j, Nz - 1, 16, Nx, Ny)] = f[idx_D3Q19(i, j, Nz - 1, 15, Nx, Ny)] + density * U;
        }
    }

    // Update macroscopic quantities
    #pragma omp parallel for schedule(static) collapse(3)
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            for (int z = 0; z < Nz; z++){
                rho[idx_D3Q19(i, j, z, Nx, Ny)] = 0;
                ux[idx_D3Q19(i, j, z, Nx, Ny)] = 0;
                uy[idx_D3Q19(i, j, z, Nx, Ny)] = 0;
                uz[idx_D3Q19(i, j, z, Nx, Ny)] = 0;
                for (int k = 0; k < 19; k++){
                    rho[idx_D3Q19(i, j, z, Nx, Ny)] += f[idx_D3Q19(i, j, z, k, Nx, Ny)];
                    ux[idx_D3Q19(i, j, z, Nx, Ny)] += c_D3Q19[k][0] * f[idx_D3Q19(i, j, z, k, Nx, Ny)];
                    uy[idx_D3Q19(i, j, z, Nx, Ny)] += c_D3Q19[k][1] * f[idx_D3Q19(i, j, z, k, Nx, Ny)];
                    uz[idx_D3Q19(i, j, z, Nx, Ny)] += c_D3Q19[k][2] * f[idx_D3Q19(i, j, z, k, Nx, Ny)];
                }
                ux[idx_D3Q19(i, j, z, Nx, Ny)] /= rho[idx_D3Q19(i, j, z, Nx, Ny)];
                uy[idx_D3Q19(i, j, z, Nx, Ny)] /= rho[idx_D3Q19(i, j, z, Nx, Ny)];
                uz[idx_D3Q19(i, j, z, Nx, Ny)] /= rho[idx_D3Q19(i, j, z, Nx, Ny)];
            }
        }
    }
};

double computeMass_3D(const std::vector<double>& rho) {
    double totalMass = 0.0;
    #pragma omp parallel for reduction(+:totalMass)
    for (int z = 0; z < Nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                int idx = (z * Ny + y) * Nx + x;
                totalMass += rho[idx];
            }
        }
    }
    return totalMass;
}
