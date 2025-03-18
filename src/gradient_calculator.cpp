#include "gradient_calculator.hpp"
#include "weight.hpp"

GradientCalculator::GradientCalculator(double& re){
    this->re  = re;
}

void GradientCalculator::calc(std::vector<Particle>& particles){
#pragma omp parallel for
    for (auto& pi : particles) {
        if(pi.type != ParticleType::Elastic)
            pi.gradient = Eigen::Vector3d::Zero();
        else{
            Eigen::Vector3d grad = Eigen::Vector3d::Zero();
            for (auto& ne : pi.neighbors){
                if (particles[ne.id].type == ParticleType::Ghost)
                    continue;
                if (particles[ne.id].type != ParticleType::Ghost){
                    Particle& pj = particles[ne.id];
                    double dij = (pi.divergence + pj.divergence)/2 ;
                    grad += dij * (pj.position - pi.position) * weight(ne.initialdistance, re) / (ne.distance * ne.initialdistance);
                }
            }
            pi.gradient = grad;
        }
    }
}

