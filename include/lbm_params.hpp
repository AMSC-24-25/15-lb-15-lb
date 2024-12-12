
constexpr int Nx = 16;
constexpr int Ny = 16;

/*
         c6   c2   c5
          \   |   /
           \  |  /
     c3 ---- c0 ---- c1
           / | \
          /  |  \
        c7   c4   c8
*/
constexpr int c[9][2] = {
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
constexpr double w[9] = {
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