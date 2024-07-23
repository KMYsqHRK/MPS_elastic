#pragma once

#include "Eigen/Dense"
#include "common.hpp"
#include <vector>

/**
 * @brief Enum class for particle type
 */
enum class ParticleType {
    Fixed,     ///< Ghost particle (outside of the domain, not used for calculation)
    Elastic,     ///< Fluid particle
    Ghost,      ///< Wall particle
};

/**
 * @brief Class for neighbor particle
 */
class Neighbor {
private:
public:
    int id;          ///< index of the neighbor particle
    double distance; ///< distance between the particle and the neighbor particle
    double initialdistance;
    double localangle = 0;
    Eigen::Vector3d strain      = Eigen::Vector3d::Zero();
    Eigen::Vector3d sheerstrain = Eigen::Vector3d::Zero();

    Neighbor(int id, double distance, double initialdist) {
        this->id       = id;
        this->distance = distance;
        this->initialdistance = initialdist;
    }
};

/**
 * @brief Class for particle in MPS method
 */
class Particle {
private:
public:
    int id;            ///< index of the particle
    ParticleType type; ///< type of the particle

    Eigen::Vector3d position;
    Eigen::Vector3d initialposition;                               ///< position of the particle
    Eigen::Vector3d velocity;                               ///< velocity of the particle
    Eigen::Vector3d acceleration = Eigen::Vector3d::Zero(); ///< acceleration of the particle
    double angle;
    double angularvelocity;
    double angurlaracceleration  = 0;
    double sourceTerm            = 0;                   ///< source term of the particle
    double divergence            = 0;
    Eigen::Vector3d gradient     = Eigen::Vector3d::Zero();
    Eigen::Vector3d Laplacian    = Eigen::Vector3d::Zero();
    double torque                = 0;

    std::vector<Neighbor> neighbors; ///< neighbors of the particle

    /**
     * @brief constructor
     * @param id  index of the particle
     * @param type  type of the particle
     * @param pos  position of the particle
     * @param vel 	velocity of the particle
     */
    Particle(int id, ParticleType type, Eigen::Vector3d pos, Eigen::Vector3d vel, double ang, double angularvel);

    /**
     * @brief calculate inverse of density
     * @param density density of the particle
     * @return inverse of density
     */
    double inverseDensity(double& density) const;
};
