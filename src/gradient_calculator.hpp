#include "particle.hpp"
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

class GradientCalculator {
public:
    GradientCalculator() = default;
    GradientCalculator(double& re);
    void calc(std::vector<Particle>& particles);

private:
    double re;
};