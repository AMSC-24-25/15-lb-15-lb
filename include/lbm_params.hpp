
/*  Used to calculate the index of D2Q9 index
    eg: f[idx_D2Q9(x, y, i, Nx)] = 1.0; 
    (x, y) means the coordinate of grid points
    i means the velocity direction at (x, y)
    Nx means the number of grid points in the x-direction
*/
inline size_t idx_D2Q9(int x, int y, int i, int Nx){
    return static_cast<size_t>((y * Nx + x) * 9 + i);
}

/*  eg: rho[idx_D2Q9(x, y, Nx)] = 1.0;
*/
inline size_t idx_D2Q9(int x, int y, int Nx){
    return static_cast<size_t>(y * Nx + x);
}

constexpr int Nx = 100;
constexpr int Ny = 100;
constexpr double tau = 0.8;
constexpr double U = 0.6;

/*
         c6   c2   c5
          \   |   /
           \  |  /
     c3 ---- c0 ---- c1
           / | \
          /  |  \
        c7   c4   c8
*/
constexpr int c_D2Q9[9][2] = {
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

//The index is the same as c
constexpr double w_D2Q9[9] = {
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