#include "divergence_calculator.hpp"
#include "weight.hpp"

DivergenceCalculator::DivergenceCalculator(int& dim, double& re){
    this->dim = dim;
    this->re  = re;
}

void DivergenceCalculator::calc(std::vector<Particle>& particles, double& n0){
#pragma omp parallel for
    for (auto& pi : particles) {
        if(pi.type == ParticleType::Ghost)
            pi.divergence = 0;
        else{
            double div = 0;
            for (auto& ne : pi.neighbors){
                if (particles[ne.id].type == ParticleType::Ghost)
                    continue;
                if (particles[ne.id].type != ParticleType::Ghost){
                    double a = ne.strain.dot(particles[ne.id].position - pi.position) * weight(ne.distance, re);
                    double b = ne.initialdistance * ne.distance;
                    div += a/b;
            }
            pi.divergence = div * this->dim / n0;
        }
    }
}
}