#ifndef LBM_PARAMS_HPP
#define LBM_PARAMS_HPP

namespace D3Q19{
    const int Q = 19;  // D3Q19速度方向数
    const int Nx = 100; 
    const int Ny = 100; // 网格尺寸 (y方向)
    const int Nz = 100; // 网格尺寸 (z方向)
    const double tau = 0.6; // 松弛时间
    const double U = 0.1;   // 顶部速度 (边界条件)
}

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

inline int idx_D3Q19(int x, int y, int z, int k, int Nx, int Ny, int Nz) {
    return ((z * D3Q19::Ny + y) * D3Q19::Nx + x) * D3Q19::Q + k;
}

#endif
