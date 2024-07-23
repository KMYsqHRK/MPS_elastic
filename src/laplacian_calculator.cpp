#include "laplacian_calculator.hpp"
#include "weight.hpp"

LaplacianCalculator::LaplacianCalculator(double& re){
    this->re  = re;
}

void LaplacianCalculator::calc(std::vector<Particle>& particles){
#pragma omp parallel for
    for (auto& pi : particles) {
        if(pi.type != ParticleType::Elastic){
            pi.Laplacian = Eigen::Vector3d::Zero();
        }
        else{
            Eigen::Vector3d laplace = Eigen::Vector3d::Zero();
            for (auto& ne : pi.neighbors){
                if (particles[ne.id].type == ParticleType::Ghost)
                    continue;
                if (particles[ne.id].type != ParticleType::Ghost){
                    double w = weight(ne.initialdistance, re) ;
                    laplace += w * ne.strain / (ne.initialdistance * ne.initialdistance);
                }
            }
            pi.Laplacian = laplace;
        }
    }
}
