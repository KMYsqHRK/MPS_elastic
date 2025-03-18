#include "particle.hpp"
Particle::Particle(int id, ParticleType type, Eigen::Vector3d pos, Eigen::Vector3d vel, double ang, double angularvel) {
    this->id       = id;
    this->type     = type;
    this->position = pos;
    this->initialposition = pos;
    this->angle = ang;
    this->angularvelocity = angularvel;
    this->velocity = vel;
}
double Particle::inverseDensity(double& density) const {
    switch (type) {
    case ParticleType::Ghost:
        return std::numeric_limits<double>::infinity();

    case ParticleType::Fixed:
        return 1 / density;

    case ParticleType::Elastic:
        return 1 / density;
    default:
        return 0;
    }
}
