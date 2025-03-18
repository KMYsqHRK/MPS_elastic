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
            double tor = 0;//トルクの初期化
            for (auto& ne : pi.neighbors) {
                if (particles[ne.id].type == ParticleType::Ghost)
                    continue;
                if (particles[ne.id].type != ParticleType::Ghost) {
                    Particle& pj = particles[ne.id];
                    Eigen::Vector3d v;
                    Eigen::Affine3d rot;
                    rot = Eigen::AngleAxisd(ne.localangle, Eigen::Vector3d(0, 0, 1));
                    v = (rot * (pj.initialposition - pi.initialposition)).cross(ne.sheerstrain);//ここの外積、怪しいな.....
                    double a = dim * lame2 * l0 * l0 * weight(ne.initialdistance, re);//ここは完全に間違ってた！！←この値が小さすぎる * weight(ne.initialdistance, re) をのぞいた。
                    double b = n0 * ne.initialdistance * ne.initialdistance;
                    tor += (v.z() * a / b);//あっていると思う。
                }
            }
            pi.torque = tor;
        }
    }
}