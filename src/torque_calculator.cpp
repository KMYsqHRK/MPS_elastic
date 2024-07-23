#include "torque_calculator.hpp"
#include "weight.hpp"

TorqueCalculator::TorqueCalculator(int& dim, double& re, double& l0, double& lame2){
    this->dim   = dim  ;
    this->re    = re   ;
    this->l0    = l0   ;
    this->lame2 = lame2;
}

void TorqueCalculator::calc(std::vector<Particle>& particles, double& n0){
#pragma omp parallel for
    for (auto& pi : particles) {
        if (pi.type == ParticleType::Ghost){
            pi.torque = 0;
        }
        else{
            double tor = 0;
            for (auto& ne : pi.neighbors) {
                if (particles[ne.id].type == ParticleType::Ghost)
                    continue;
                if (particles[ne.id].type != ParticleType::Ghost) {
                    Particle& pj = particles[ne.id];
                    Eigen::Vector3d v;
                    Eigen::Affine3d rot;
                    rot = Eigen::AngleAxisd(ne.localangle, Eigen::Vector3d(0, 0, 1));
                    v = -1 * (rot * (pj.initialposition - pi.initialposition)).cross(ne.sheerstrain);
                    double a = this->dim * this->lame2 * weight(ne.initialdistance, this->re) * this->l0 * this->l0 * this->l0;
                    double b = n0 * ne.initialdistance * ne.initialdistance;
                    tor -= v.z() * a / b;
                }
            }
            pi.torque = tor;
        }
    }
}