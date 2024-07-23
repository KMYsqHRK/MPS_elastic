#include "particle.hpp"
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

class TorqueCalculator {
public:
    TorqueCalculator() = default;
    TorqueCalculator(int& dim, double& re, double& l0, double& lame2);
    void calc(std::vector<Particle>& particles, double& n0);

private:
    int dim;
    double re;
    double l0;
    double lame2;
};
