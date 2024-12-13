#include <vector>

void D2Q9_serial(std::vector<double>& f, std::vector<double>& rho, std::vector<double>& ux, std::vector<double>& uy, std::vector<double>& f_eq);
void D2Q9_parallel(std::vector<double>& f, std::vector<double>& rho, std::vector<double>& ux, std::vector<double>& uy, std::vector<double>& f_eq);
