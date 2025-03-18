#include "particle.hpp"
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

class Hamiltonian {
public:
    Hamiltonian() = default;
    Hamiltonian(int& dim, double& re, double& lame1, double& lame2, double& rho);
    void calc(std::vector<Particle>& particles);
    void update_acc(std::vector<Particle>& particles);
private:
    int dim;
    double re;
    double lame1;
    double lame2; 
    double rho;  
};