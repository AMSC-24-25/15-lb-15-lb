// lbm_params.hpp
#ifndef LBM_PARAMS_HPP
#define LBM_PARAMS_HPP

#include <cmath>

// 域大小
const int Nx = 200; // x方向的网格点数
const int Ny = 100; // y方向的网格点数

// LBM 参数
const double U = 0.1;                // 盖子的速度
const double Re = 100.0;             // Reynolds 数
const double rho0 = 1.0;             // 初始密度
const double obstacle_radius = 10.0; // 圆形障碍物的半径
const int obstacle_center_x = Nx / 4; // 障碍物中心的 x 坐标
const int obstacle_center_y = Ny / 2; // 障碍物中心的 y 坐标

// 计算粘度和松弛时间 tau
const double nu = U * obstacle_radius / Re;    // 动力粘度
const double tau = 3.0 * nu + 0.5;             // 松弛时间

// D2Q9 格子速度
const int c_D2Q9[9][2] = {
    {0, 0},  // 0: 静止
    {1, 0},  // 1: 东
    {0, 1},  // 2: 北
    {-1, 0}, // 3: 西
    {0, -1}, // 4: 南
    {1, 1},  // 5: 东北
    {-1, 1}, // 6: 西北
    {-1, -1},// 7: 西南
    {1, -1}  // 8: 东南
};

// D2Q9 权重
const double w_D2Q9[9] = {
    4.0 / 9.0,     // 0
    1.0 / 9.0,     // 1
    1.0 / 9.0,     // 2
    1.0 / 9.0,     // 3
    1.0 / 9.0,     // 4
    1.0 / 36.0,    // 5
    1.0 / 36.0,    // 6
    1.0 / 36.0,    // 7
    1.0 / 36.0     // 8
};

// 计算 D2Q9 索引
inline int idx_D2Q9(int i, int j, int k, int Nx) {
    return (j * Nx + i) * 9 + k;
}

inline int idx_D2Q9(int i, int j, int Nx) {
    return j * Nx + i;
}

#endif // LBM_PARAMS_HPP