#include "particle.hpp"
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

class DivergenceCalculator {
public:
    DivergenceCalculator() = default;
    DivergenceCalculator(int& dim, double& re);
    void calc(std::vector<Particle>& particles, double& n0);

private:
    int dim;
    double re;
};