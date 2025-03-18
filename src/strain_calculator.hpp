#include "particle.hpp"
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Geometry>

class StrainCalculator {
public:
    StrainCalculator() = default;
    StrainCalculator(double& re);
    void calc(std::vector<Particle>& particles);

private:
    double re;
};
