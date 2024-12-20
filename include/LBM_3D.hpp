#include <string>
#include <vector>

class LBM_3D {
  public:
    // Constructor
    LBM_3D(int Nx, int Ny, int Nz, double U, double rho0);

    // Public function
    void initialize();
    void D3Q19_parallel_iterate(int step);
    void save_to_CSV(const std::string &filename) const;

  private:
    int Nx, Ny, Nz;
    double U, rho0;

    constexpr static int Q = 19;
    constexpr static double tau = 0.6;
    constexpr static double omega = 1.0 / tau;

    constexpr static int c_D3Q19[19][3] = {
    {0, 0, 0}, {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1},
    {1, 1, 0}, {-1, -1, 0}, {1, -1, 0}, {-1, 1, 0}, {1, 0, 1}, {-1, 0, -1},
    {1, 0, -1}, {-1, 0, 1}, {0, 1, 1}, {0, -1, -1}, {0, 1, -1}, {0, -1, 1}
    };


    constexpr static double w_D3Q19[19] = {
        1.0 / 3.0,   // 0
        1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
    };

    // Velocity distribution function
    std::vector<double> f;
    std::vector<double> f_temp;

    // Density
    std::vector<double> rho;

    // Velocity in x-direction
    std::vector<double> ux;

    // Velocity in y-direction
    std::vector<double> uy;

    //Velocity in z-direction
    std::vector<double> uz;

    inline size_t idx_D3Q19(int i, int j, int z, int k, int Nx, int Ny) {
        return static_cast<size_t>(((z * Ny + j) * Nx + i) * 19 + k);
    }

    inline size_t idx_D3Q19(int i, int j, int z, int Nx, int Ny) {
        return static_cast<size_t> ((z * Ny + j) * Nx + i);
    }
};