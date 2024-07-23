#include "strain_calculator.hpp"
#include "weight.hpp"

StrainCalculator::StrainCalculator(double& re) {
    this->re = re;
}

//calculate localangle, strain, sheerstrain
void StrainCalculator::calc(std::vector<Particle>& particles) {
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Ghost)
            continue;
        for (auto& ne : pi.neighbors) {
            if (particles[ne.id].type == ParticleType::Ghost)
                continue;
            if (particles[ne.id].type != ParticleType::Ghost) {
                Particle& pj = particles[ne.id];
                double Localangle = (pi.angle + pj.angle)/2;
                Eigen::Affine3d rot;
                rot = Eigen::AngleAxisd(Localangle, Eigen::Vector3d(0, 0, 1));
                ne.strain  =  
                (pj.position - pi.position) - rot *(pj.initialposition - pi.initialposition);
                ne.localangle = Localangle;
                double a = ne.strain.dot(pj.position - pi.position);
                Eigen::Vector3d normalstrain = 
                (pj.position - pi.position) * a / (ne.distance * ne.distance);
                ne.sheerstrain = ne.strain - normalstrain;
            }
        }
    }
}