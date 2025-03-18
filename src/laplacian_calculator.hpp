#include "particle.hpp"
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

class LaplacianCalculator {
public:
    LaplacianCalculator() = default;
    LaplacianCalculator(double& re);
    void calc(std::vector<Particle>& particles);

private:
    double re;
};