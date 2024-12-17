#include <string>
#include <vector>

class LBM_2D {
  public:
    // Constructor
    LBM_2D(int Nx, int Ny, int Re, double U, double rho0);

    // Public functions
    void initialize();
    void D2Q9_serial_iterate(int step);
    void D2Q9_parallel_iterate(int step);
    void save_to_CSV(const std::string &filename) const;

  private:
    int Nx, Ny, Re;
    double U, rho0;

    constexpr static double c_s = 0.5773502692;
    constexpr static double c_s_square = c_s * c_s;
    constexpr static double nu = 0.0667;
    constexpr static double tau = nu / (c_s * c_s) + 0.5;

    /*
     c6   c2   c5
      \   |   /
       \  |  /
 c3 ---- c0 ---- c1
       / | \
      /  |  \
    c7   c4   c8
    */
    constexpr static int c_D2Q9[9][2] = {
        {0, 0},
        {1, 0},
        {0, 1},
        {-1, 0},
        {0, -1},
        {1, 1},
        {-1, 1},
        {-1, -1},
        {1, -1}
    };

    // The index is the same as c
    constexpr static double w_D2Q9[9] = {
        4.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0
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

    /*  Used to calculate the index of D2Q9 index
    eg: f[idx_D2Q9(x, y, i, Nx)] = 1.0;
    (x, y) means the coordinate of grid points
    i means the velocity direction at (x, y)
    Nx means the number of grid points in the x-direction
    */
    inline size_t idx_D2Q9(int x, int y, int i, int Nx) {
        return static_cast<size_t>((y * Nx + x) * 9 + i);
    }

    // eg: rho[idx_D2Q9(x, y, Nx)] = 1.0;
    inline size_t idx_D2Q9(int x, int y, int Nx) {
        return static_cast<size_t>(y * Nx + x);
    }
};