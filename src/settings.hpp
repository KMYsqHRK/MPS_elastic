#pragma once

#include "common.hpp"
#include "domain.hpp"

#include <Eigen/Dense>
#include <filesystem>

/**
 * @brief Struct for settings of calculation
 *
 * @details This struct contains the settings for the simulation. It is used to
 * load the settings from the input file.
 */
struct Settings {
    // computational condition
    int dim{};                 ///< Dimension of the simulation
    double particleDistance{}; ///< Initial distance between particles
    double dt{};               ///< Time step
    double endTime{};          ///< End time of the simulation
    double outputPeriod{};     ///< Output period of the simulation
    double cflCondition{};     ///< CFL condition
    int numPhysicalCores{};    ///< Number of cores to calculate

    Domain domain{}; ///< domain of the simulation

    double elasticdensity{};
    double elasticratio{}; 
    double poissonratio{}; 
    double lame1{};
    double lame2{};

    // gravity
    Eigen::Vector3d gravity; ///< Gravity

    // collision
    double collisionDistance{};        ///< Distance for collision detection
    double coefficientOfRestitution{}; ///< Coefficient of restitution

    // effective radius
    double re_forDivergence{}; ///< Effective radius for number density
    double re_forGradient{};      ///< Effective radius for gradient
    double re_forLaplacian{};     ///< Effective radius for Laplacian
    double reMax{};               ///< Maximum of effective radius

    // i/o
    std::filesystem::path profPath; ///< Path for input particle file
};
