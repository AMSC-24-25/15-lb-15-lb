#include <vector>

void D2Q9_serial(std::vector<double>& f, std::vector<double>& rho, std::vector<double>& ux, std::vector<double>& uy, std::vector<double>& f_eq);
void D2Q9_parallel(std::vector<double>& f, std::vector<double>& rho, std::vector<double>& ux, std::vector<double>& uy, std::vector<double>& f_eq);

void initializeData(std::vector<double>& f, std::vector<double>& rho, std::vector<double>& ux, std::vector<double>& uy, std::vector<double>& uz);

void D3Q19_parallel(std::vector<double>& f, std::vector<double>& rho, std::vector<double>& ux, std::vector<double>& uy, std::vector<double>& uz, std::vector<double>& f_eq);

double computeMass_2D(const std::vector<double>& rho);

double computeMass_3D(const std::vector<double>& rho);



void saveToCSV_D2Q9(const std::string& filename, int Nx, int Ny, 
               const std::vector<double>& rho, 
               const std::vector<double>& ux, 
               const std::vector<double>& uy);

void saveToCSV_D3Q19(const std::string& filename, const std::vector<double>& rho,
               const std::vector<double>& ux, const std::vector<double>& uy,
               const std::vector<double>& uz);
