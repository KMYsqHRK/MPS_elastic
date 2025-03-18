#pragma once

#include "bucket.hpp"
#include "common.hpp"
#include "domain.hpp"
#include "input.hpp"
#include "neighbor_searcher.hpp"
#include "divergence_calculator.hpp"
#include "torque_calculator.hpp"
#include "strain_calculator.hpp"
#include "gradient_calculator.hpp"
#include "laplacian_calculator.hpp"
#include "refvalues.hpp"
#include "settings.hpp"
#include "Hamiltonian.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <vector>

/***
 * @brief MPS simulation class
 *
 * @details Executes the MPS simulation. This class does not handle the
 * simulation process itself, but only the calculation of the MPS method.
 */
class MPS {
public:
    Settings settings;                   ///< Settings for the simulation
    RefValues refValuesForDivergence; ///< Reference values for the simulation (\f$n^0\f$, \f$\lambda^0\f$)
    RefValues refValuesForLaplacian;     ///< Reference values for the simulation (\f$n^0\f$, \f$\lambda^0\f$)
    RefValues refValuesForGradient;      ///< Reference values for the simulation (\f$n^0\f$, \f$\lambda^0\f$)
    std::vector<Particle> particles;     ///< Particles in the simulation
    Domain domain{};                     ///< Domain of the simulation

    double courant{}; ///< Maximum courant number among all particles

    MPS() = default;

    MPS(const Input& input);

    void stepForward(int N);

private:
    DivergenceCalculator divergencecalculator;
    TorqueCalculator torquecalculator;
    StrainCalculator straincalculator;
    GradientCalculator gradientcalculator;
    LaplacianCalculator laplaciancalculator;
    NeighborSearcher neighborSearcher;
    Hamiltonian hamiltonian;

    void caldivergenceGradient(const double& re);
    void moveParticle();
    void collision();
    void calgravity();
    void rotateParticle();
    void calcstrainLaplacian(const double& re);
    void calacceleration();
    void calangularacceleration();
    void upgradeNeighbors();
    /**
     * @brief calculate Courant number
     */
    void calCourant();
};
