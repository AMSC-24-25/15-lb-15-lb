// main.cpp
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <sys/time.h>
#include "lbm_params.hpp"
#include "function.hpp"

// 初始化函数：设置固体节点和初始分布函数
void initialize(std::vector<double>& f, 
                std::vector<double>& rho, 
                std::vector<double>& ux, 
                std::vector<double>& uy, 
                std::vector<int>& is_solid) {
    // 初始化固体节点（圆形障碍物）
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            double dx = i - obstacle_center_x;
            double dy = j - obstacle_center_y;
            if (std::sqrt(dx * dx + dy * dy) <= obstacle_radius) {
                is_solid[j * Nx + i] = 1;
            }
        }
    }

    // 初始化分布函数为平衡分布
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            if (is_solid[j * Nx + i]) {
                // 固体节点：设置为静止分布
                for (int k = 0; k < 9; k++) {
                    f[idx_D2Q9(i, j, k, Nx)] = w_D2Q9[k] * rho0;
                }
                continue;
            }

            // 设置顶边界（盖子）的速度
            if (j == Ny - 1) {
                ux[j * Nx + i] = U;
                uy[j * Nx + i] = 0.0;
            }

            // 计算平衡分布函数
            for (int k = 0; k < 9; k++) {
                double cu = c_D2Q9[k][0] * ux[j * Nx + i] + c_D2Q9[k][1] * uy[j * Nx + i];
                double u_sq = ux[j * Nx + i] * ux[j * Nx + i] + uy[j * Nx + i] * uy[j * Nx + i];
                f[idx_D2Q9(i, j, k, Nx)] = w_D2Q9[k] * rho0 * 
                    (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq);
            }
        }
    }
}

int main(){
    // 初始化分布函数、密度、速度等
    std::vector<double> f(Nx * Ny * 9, 0.0);
    std::vector<double> rho(Nx * Ny, rho0);
    std::vector<double> ux(Nx * Ny, 0.0);
    std::vector<double> uy(Nx * Ny, 0.0);
    std::vector<double> f_eq(Nx * Ny * 9, 0.0);

    // 初始化固体标志数组
    std::vector<int> is_solid(Nx * Ny, 0);
    initialize(f, rho, ux, uy, is_solid);

    // 力的累积变量
    double total_drag = 0.0;
    double total_lift = 0.0;

    // 计时变量
    struct timeval t1, t2;
    double etime;

    // 开始模拟循环
    gettimeofday(&t1, NULL);
    int max_iters = 10000;
    for (int iter = 0; iter < max_iters; iter++) {
        double iter_drag = 0.0;
        double iter_lift = 0.0;

        // 执行并行的 LBM 步骤
        D2Q9_parallel(f, rho, ux, uy, f_eq, is_solid, iter_drag, iter_lift);
        total_drag += iter_drag;
        total_lift += iter_lift;

        // 施加顶边界的速度（Zou/He 边界条件）
        #pragma omp parallel for
        for(int i = 0; i < Nx; i++) {
            int j = Ny - 1;
            // 重新计算 rho
            double rho_lid = f[idx_D2Q9(i, j, 0, Nx)] + f[idx_D2Q9(i, j, 1, Nx)] + f[idx_D2Q9(i, j, 3, Nx)] 
                            + 2.0 * (f[idx_D2Q9(i, j, 2, Nx)] + f[idx_D2Q9(i, j, 5, Nx)] + f[idx_D2Q9(i, j, 8, Nx)]);
            rho[j * Nx + i] = rho_lid;

            // 设定速度
            ux[j * Nx + i] = U;
            uy[j * Nx + i] = 0.0;

            // 更新分布函数（按照 Zou/He 边界条件）
            // f4 = f2
            f[idx_D2Q9(i, j, 4, Nx)] = f[idx_D2Q9(i, j, 2, Nx)];

            // f7 = f5 + 0.5*(f3 - f1) + (1.0 / 12.0)*rho*U
            f[idx_D2Q9(i, j, 7, Nx)] = f[idx_D2Q9(i, j, 5, Nx)] + 
                                       0.5 * (f[idx_D2Q9(i, j, 3, Nx)] - f[idx_D2Q9(i, j, 1, Nx)]) + 
                                       (1.0 / 12.0) * rho[j * Nx + i] * ux[j * Nx + i];

            // f6 = f8 + 0.5*(f1 - f3)
            f[idx_D2Q9(i, j, 6, Nx)] = f[idx_D2Q9(i, j, 8, Nx)] + 
                                       0.5 * (f[idx_D2Q9(i, j, 1, Nx)] - f[idx_D2Q9(i, j, 3, Nx)]);
        }

        // 质量守恒检查
        if ((iter + 1) % 1000 == 0) {
            double total_mass = 0.0;
            #pragma omp parallel for reduction(+:total_mass)
            for(int idx = 0; idx < rho.size(); idx++) {
                total_mass += rho[idx];
            }
            double expected_mass = rho0 * Nx * Ny;
            printf("Iteration %d: Drag = %lf, Lift = %lf, Total Mass = %lf (Expected: %lf)\n", 
                   iter + 1, iter_drag, iter_lift, total_mass, expected_mass);
        }
    }
    gettimeofday(&t2, NULL);

    // 计算并打印总耗时
    etime = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
    etime /= 1000.0;
    printf("Parallel simulation completed in %f seconds.\n", etime);
    printf("Total Drag: %lf\n", total_drag);
    printf("Total Lift: %lf\n", total_lift);

    // 验证：中心点的速度和密度
    int center_i = Nx / 2;
    int center_j = Ny / 2;
    double center_rho = rho[idx_D2Q9(center_i, center_j, Nx)];
    double center_ux_val = ux[idx_D2Q9(center_i, center_j, Nx)];
    double center_uy_val = uy[idx_D2Q9(center_i, center_j, Nx)];
    printf("Center Density: %lf\n", center_rho);
    printf("Center Velocity Ux: %lf\n", center_ux_val);
    printf("Center Velocity Uy: %lf\n", center_uy_val);

    // 保存结果到 CSV 文件
    saveToCSV("lbm_results.csv", Nx, Ny, rho, ux, uy);
    printf("Results saved to lbm_results.csv\n");

    return 0;
}