#include "mps.hpp"

#include "particle.hpp"
#include "weight.hpp"

#include <queue>

using std::cerr;
using std::endl;

/**
 * @brief Initialize settings 
 * @details Initialize settings variable and prepare some values for calculating MPS method.
 * constructor of MPS method. 
 */
MPS::MPS(const Input& input) {
    this->settings             = input.settings;
    this->domain               = input.settings.domain;
    this->particles            = input.particles;
    refValuesForDivergence     = RefValues(settings.dim, settings.particleDistance, settings.re_forDivergence);
    refValuesForGradient       = RefValues(settings.dim, settings.particleDistance, settings.re_forGradient);
    refValuesForLaplacian      = RefValues(settings.dim, settings.particleDistance, settings.re_forLaplacian);
    this->divergencecalculator = DivergenceCalculator(settings.dim, settings.re_forDivergence);
    this->torquecalculator     = TorqueCalculator(settings.dim, settings.re_forDivergence, settings.particleDistance, settings.lame2);
    this->straincalculator     = StrainCalculator(settings.re_forDivergence);
    this->gradientcalculator   = GradientCalculator(settings.re_forGradient);
    this->laplaciancalculator  = LaplacianCalculator(settings.re_forLaplacian);
    this->neighborSearcher     = NeighborSearcher(settings.reMax, settings.domain, particles.size());
    this->hamiltonian          = Hamiltonian(settings.dim, settings.particleDistance, settings.lame1, settings.lame2, settings.elasticdensity);
}

/**
 * @brief Process simulation
 * @details First of all, gravity and forces move particles into temporary positions. 
 * Next, the position of the colliding particle is adjusted, the particle number density is calculated and the particle boundary conditions are redefined. 
 * Next, the minimum pressure is found, the pressure gradient is calculated, and the position of the particle is determined according to the pressure gradient.
 */
void MPS::stepForward(int N) {
    /// only at first step, set neighbors 
    if (N == 0) {
        neighborSearcher.setNeighbors(particles); //検証済み。
        
    }
    ///隣接粒子情報を更新する
    upgradeNeighbors(); //検証済み。

    ///ここ以下は支配方程式次第で必要な演算が異なる。(上がMPS、下がHamiltonian)
    if (1 == 1) {
        straincalculator.calc(particles);

        divergencecalculator.calc(particles, refValuesForDivergence.n0);

        gradientcalculator.calc(particles);

        laplaciancalculator.calc(particles);

        torquecalculator.calc(particles, refValuesForDivergence.n0);

        calgravity();

        calangularacceleration(); 

        rotateParticle();

        calacceleration();

        moveParticle();

        ///collision();

        calCourant();
    }
    if(1 == 0){
        calgravity();

        hamiltonian.calc(particles);

        hamiltonian.update_acc(particles);

        moveParticle();

        calCourant();
    }
}

/**
 * @brief Calculate gravitational acceleration and add to acceleration
 */
void MPS::calgravity() {
#pragma omp parallel for
    for (auto& p : particles) {
        if (p.type == ParticleType::Elastic) {
            p.acceleration = settings.gravity;//加速度の初期化
        }
    }
}

void MPS::calacceleration() {
#pragma omp parallel for
    for (auto& p : particles){
        if (p.type == ParticleType::Elastic) {
            p.acceleration += (settings.lame1 * p.gradient / refValuesForGradient.n0 + 2 * settings.lame2 * p.Laplacian / refValuesForLaplacian.n0) * settings.dim / settings.elasticdensity;
        } 
    }
}

void MPS::calangularacceleration() {
#pragma omp parallel for
    for (auto& p : particles){
        p.angurlaracceleration = 0;//回転加速度の初期化
        if (p.type != ParticleType::Ghost) {
            p.angurlaracceleration += p.torque * 6 /(settings.elasticdensity * settings.particleDistance * settings.particleDistance * settings.particleDistance * settings.particleDistance);
        } 
    }
}

void MPS::rotateParticle() {
#pragma omp parallel for
    for (auto& p : particles) {
        if (p.type != ParticleType::Ghost) {
            p.angularvelocity += p.angurlaracceleration * settings.dt;
            p.angle += p.angularvelocity * settings.dt; //各速度を計算してから角度を計算する。(安定しやすくなるかもしれない?)
            //p.angularvelocity = p.angurlaracceleration * settings.dt;
        }
        p.angurlaracceleration = 0; //回転加速度の初期化
    }
}

/**
 * @brief Move particles using temporal acceleration
 */
void MPS::moveParticle() {
#pragma omp parallel for
    for (auto& p : particles) {
        if (p.type == ParticleType::Elastic) {
            p.velocity += p.acceleration * settings.dt;
            p.position += p.velocity * settings.dt;//速度を計算してから、位置を計算する。
            //p.velocity += p.acceleration * settings.dt;
        }
    }
}

/**
 * @brief If particles get closer than collision distance, perform collision processing
 */
void MPS::collision() {
    for (auto& pi : particles) {
        if (pi.type != ParticleType::Elastic)
            continue;

        for (auto& neighbor : pi.neighbors) {
            Particle& pj = particles[neighbor.id];
            if (pj.type == ParticleType::Elastic && pj.id >= pi.id)
                continue;

            if (neighbor.distance < settings.collisionDistance) {

                double invMi = pi.inverseDensity(settings.elasticdensity);
                double invMj = pj.inverseDensity(settings.elasticdensity);
                double mass  = 1.0 / (invMi + invMj);

                Eigen::Vector3d normal = (pj.position - pi.position).normalized();
                double relVel          = (pj.velocity - pi.velocity).dot(normal);
                double impulse         = 0.0;
                if (relVel < 0.0)
                    impulse = -(1 + settings.coefficientOfRestitution) * relVel * mass;
                pi.velocity -= impulse * invMi * normal;
                pj.velocity += impulse * invMj * normal;

                double depth           = settings.collisionDistance - neighbor.distance;
                double positionImpulse = depth * mass;
                pi.position -= positionImpulse * invMi * normal;
                pj.position += positionImpulse * invMj * normal;

                // cerr << "WARNING: Collision between particles " << pi.id << " and " << pj.id << " occurred."
                // << endl;
            }
        }
    }
}



/**
 * @brief Check for errors with courant number
 */
void MPS::calCourant() {
    courant = 0.0;

    for (auto& pi : particles) {
        if (pi.type != ParticleType::Elastic)
            continue;

        double iCourant = (pi.velocity.norm() * settings.dt) / settings.particleDistance;
        courant         = std::max(courant, iCourant);
    }

    if (courant > settings.cflCondition) {
        cerr << "ERROR: Courant number is larger than CFL condition. Courant = " << courant << endl;
    }
}
//Ghostでないすべての粒子に対して隣接粒子の情報を書き換える。
void MPS::upgradeNeighbors() {
#pragma omp parallel for
    for (auto& pi : particles){
        if (pi.type == ParticleType::Ghost)
            continue;
        for (auto& ne : pi.neighbors){
            /*
            //デバッグ用(隣接粒子を取得できているのかを確認):確認済み。隣接粒子は取得できている。
            if (pi.id == 300){
                std::cout << "neighbor : " << ne.id << std::endl;
            }
            */
            Particle& pj = particles[ne.id];
            ne.distance = (pj.position - pi.position).norm();
            ne.localangle = (pi.angle + pj.angle)/2;
        }
    }
}
