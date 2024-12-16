#ifndef LBM_PARAMS_HPP
#define LBM_PARAMS_HPP

const int Q = 19;  // D3Q19速度方向数
const int Nx = 100; 
const int Ny = 100; // 网格尺寸 (y方向)
const int Nz = 100; // 网格尺寸 (z方向)
const double tau = 0.6; // 松弛时间
const double U = 0.1;   // 顶部速度 (边界条件)
const double omega = 1.0 / tau;             // 弛豫因子

const int c_D3Q19[19][3] = {
    {0, 0, 0}, {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1},
    {1, 1, 0}, {-1, -1, 0}, {1, -1, 0}, {-1, 1, 0}, {1, 0, 1}, {-1, 0, -1},
    {1, 0, -1}, {-1, 0, 1}, {0, 1, 1}, {0, -1, -1}, {0, 1, -1}, {0, -1, 1}
};


const double w_D3Q19[19] = {
    1.0 / 3.0,   // 0
    1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
};

inline size_t idx_D3Q19(int i, int j, int z, int k, int Nx, int Ny) {
    return static_cast<size_t>(((z * Ny + j) * Nx + i) * 19 + k);
}

inline size_t idx_D3Q19(int i, int j, int z, int Nx, int Ny) {
    return static_cast<size_t> ((z * Ny + j) * Nx + i);
}


#endif
