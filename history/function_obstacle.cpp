// function.cpp
#include "function.hpp"
#include "lbm_params.hpp"
#include <omp.h>
#include <cmath>
#include <fstream>
#include <limits>

// 并行 D2Q9 LBM 步骤的实现
void D2Q9_parallel(std::vector<double>& f, 
                   std::vector<double>& rho, 
                   std::vector<double>& ux, 
                   std::vector<double>& uy, 
                   std::vector<double>& f_eq, 
                   const std::vector<int>& is_solid,
                   double& drag, 
                   double& lift) {
    double drag_sum = 0.0;
    double lift_sum = 0.0;

    // 碰撞和力计算步骤
    #pragma omp parallel for collapse(2) reduction(+:drag_sum, lift_sum)
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            // 检查当前节点是否为固体节点
            if (is_solid[j * Nx + i]) {
                // 反弹回来的碰撞：反射分布函数
                for (int k = 0; k < 9; ++k) {
                    int opposite_k;
                    switch(k) {
                        case 1: opposite_k = 3; break;
                        case 2: opposite_k = 4; break;
                        case 3: opposite_k = 1; break;
                        case 4: opposite_k = 2; break;
                        case 5: opposite_k = 7; break;
                        case 6: opposite_k = 8; break;
                        case 7: opposite_k = 5; break;
                        case 8: opposite_k = 6; break;
                        default: opposite_k = 0; break; // case 0 remains 0
                    }

                    // 动量交换计算，用于力的计算
                    if (k != 0) { // 跳过静止粒子
                        double f_k = f[idx_D2Q9(i, j, k, Nx)];
                        double f_opposite = f[idx_D2Q9(i, j, opposite_k, Nx)];
                        double delta_f = f_k - f_opposite;
                        double fx = delta_f * c_D2Q9[k][0];
                        double fy = delta_f * c_D2Q9[k][1];
                        drag_sum += fx;
                        lift_sum += fy;
                    }

                    // 应用反射
                    f[idx_D2Q9(i, j, k, Nx)] = f[idx_D2Q9(i, j, opposite_k, Nx)];
                }
                continue; // 跳过固体节点的流体计算
            }

            // 计算宏观变量
            double local_rho = 0.0;
            double local_ux_temp = 0.0;
            double local_uy_temp = 0.0;

            for (int k = 0; k < 9; ++k) {
                local_rho += f[idx_D2Q9(i, j, k, Nx)];
                local_ux_temp += f[idx_D2Q9(i, j, k, Nx)] * c_D2Q9[k][0];
                local_uy_temp += f[idx_D2Q9(i, j, k, Nx)] * c_D2Q9[k][1];
            }

            // 防止除以零或负密度
            if (local_rho <= 0.0) {
                local_rho = rho0;
            }

            rho[j * Nx + i] = local_rho;
            ux[j * Nx + i] = local_ux_temp / local_rho;
            uy[j * Nx + i] = local_uy_temp / local_rho;

            // 计算平衡分布函数
            for (int k = 0; k < 9; ++k) {
                double cu = c_D2Q9[k][0] * ux[j * Nx + i] + c_D2Q9[k][1] * uy[j * Nx + i];
                double u_sq = ux[j * Nx + i] * ux[j * Nx + i] + uy[j * Nx + i] * uy[j * Nx + i];
                f_eq[idx_D2Q9(i, j, k, Nx)] = w_D2Q9[k] * rho[j * Nx + i] * 
                    (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq);
            }

            // 碰撞步骤
            for (int k = 0; k < 9; ++k) {
                f[idx_D2Q9(i, j, k, Nx)] += -(f[idx_D2Q9(i, j, k, Nx)] - f_eq[idx_D2Q9(i, j, k, Nx)]) / tau;
            }
        }
    }

    // 流动步骤：传播分布函数
    std::vector<double> f_temp(f.size(), 0.0);

    #pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            for (int k = 0; k < 9; ++k) {
                int ni = i + c_D2Q9[k][0];
                int nj = j + c_D2Q9[k][1];

                // 处理域边界：反弹回来的边界条件
                if (ni < 0 || ni >= Nx || nj < 0 || nj >= Ny) {
                    // 各向同性反射
                    int opposite_k;
                    switch(k) {
                        case 1: opposite_k = 3; break;
                        case 2: opposite_k = 4; break;
                        case 3: opposite_k = 1; break;
                        case 4: opposite_k = 2; break;
                        case 5: opposite_k = 7; break;
                        case 6: opposite_k = 8; break;
                        case 7: opposite_k = 5; break;
                        case 8: opposite_k = 6; break;
                        default: opposite_k = 0; break; // case 0 remains 0
                    }

                    f_temp[idx_D2Q9(i, j, opposite_k, Nx)] += f[idx_D2Q9(i, j, k, Nx)];
                }
                else {
                    if (is_solid[nj * Nx + ni]) {
                        // 固体边界的反射
                        int opposite_k;
                        switch(k) {
                            case 1: opposite_k = 3; break;
                            case 2: opposite_k = 4; break;
                            case 3: opposite_k = 1; break;
                            case 4: opposite_k = 2; break;
                            case 5: opposite_k = 7; break;
                            case 6: opposite_k = 8; break;
                            case 7: opposite_k = 5; break;
                            case 8: opposite_k = 6; break;
                            default: opposite_k = 0; break; // case 0 remains 0
                        }

                        f_temp[idx_D2Q9(i, j, opposite_k, Nx)] += f[idx_D2Q9(i, j, k, Nx)];
                    }
                    else {
                        // 正常传播到邻居节点
                        f_temp[idx_D2Q9(ni, nj, k, Nx)] += f[idx_D2Q9(i, j, k, Nx)];
                    }
                }
            }
        }
    }

    // 更新分布函数
    #pragma omp parallel for
    for (size_t idx = 0; idx < f.size(); ++idx) {
        f[idx] = f_temp[idx];
    }

    // 累加总阻力和升力
    drag += drag_sum;
    lift += lift_sum;
}

// 将模拟结果保存到 CSV 文件
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