#include "weight.hpp"
#include "particle.hpp"
#include "Hamiltonian.hpp"

#include <Eigen/Dense>

Hamiltonian::Hamiltonian(int& dim, double& re, double& lame1, double& lame2, double& rho){
    this->dim = dim;
    this->re  = re;
    this->lame1 = lame1;
    this->lame2 = lame2;
    this->rho = rho;
}

void Hamiltonian::calc(std::vector<Particle>& particles){
#pragma omp parallel for
    for (auto& pi : particles){
        Eigen::Matrix3d Bi = Eigen::Matrix3d::Zero();
        for (auto& ne : pi.neighbors){
            Eigen::Vector3d realtive_pos = particles[ne.id].position - pi.position;
            Eigen::Vector3d initial_realtive_pos = particles[ne.id].initialposition - pi.initialposition;
            pi.Ai = pi.Ai + initial_realtive_pos * initial_realtive_pos.transpose() * weight(ne.initialdistance, re);
            Bi = Bi + realtive_pos * initial_realtive_pos.transpose() * weight(ne.initialdistance, re);
        }
        pi.Ai(2,2) = 1;
        Bi(2,2) = 1;
        Eigen::Matrix3d Fi = Bi * pi.Ai.inverse();
        Eigen::Matrix3d Ei = 0.5 * (Fi.transpose() * Fi - Eigen::Matrix3d::Identity());
        Eigen::Matrix3d Si =  
            lame1 * Ei.trace() * Eigen::Matrix3d::Identity() + 
            2 * lame2 * Ei;
        pi.PKstress = Fi * Si;
    }
}

void Hamiltonian::update_acc(std::vector<Particle>& particles){
#pragma omp parallel for
    for (auto& pi : particles){
        if(pi.type == ParticleType::Elastic){
            Eigen::Vector3d Vt = Eigen::Vector3d::Zero();
            for (auto& ne : pi.neighbors){
                double W = weight(ne.initialdistance, re);
                Eigen::Vector3d realtive_pos = particles[ne.id].position - pi.position;
                Eigen::Vector3d initial_realtive_pos = particles[ne.id].initialposition - pi.initialposition;
                Vt += (pi.PKstress * pi.Ai.inverse() * initial_realtive_pos + particles[ne.id].PKstress * particles[ne.id].Ai.inverse() * initial_realtive_pos) * W;
            }
            pi.acceleration += Vt/rho;
        }
        pi.PKstress.setZero();
        pi.Ai.setZero();//initialize PJstress and Ai
    }
}