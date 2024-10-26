/**
 * @file      nbody.cu
 *
 * @author    Jan Hol�� \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            xholan11@fit.vutbr.cz
 *
 * @brief     PCG Assignment 1
 *
 * @version   2024
 *
 * @date      20 October   2024, 09:00 \n
 */

#include <device_launch_parameters.h>

#include "nbody.cuh"

/* Constants */
constexpr float G = 6.67384e-11f;
constexpr float COLLISION_DISTANCE = 0.01f;

// create min float constant used on GPU (because this construction cannot be used in __global__ function)
__device__ __constant__ float MIN_FLOAT = std::numeric_limits<float>::min();

/**
 * CUDA kernel to calculate gravitation velocity
 * @param p      - particles
 * @param tmpVel - temp array for velocities
 * @param N      - Number of particles
 * @param dt     - Size of the time step
 */
__global__ void calculateGravitationVelocity(Particles p, Velocities tmpVel, const unsigned N, float dt)
{
    // calculate thread position in N particles
    const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

    // if out of bounds, end thread work
    if (idx >= N)
        return;

    // params of the current thread's particle
    const float posX = p.posX[idx];
    const float posY = p.posY[idx];
    const float posZ = p.posZ[idx];
    const float weight = p.weight[idx];

    // place where the new params will be stored
    float newVelX{};
    float newVelY{};
    float newVelZ{};

    const float G_dt = G * dt;
    // loop through all other particles and compute the differences
    for (unsigned j = 0u; j < N; ++j)
    {
        const float dx = p.posX[j] - posX;
        const float dy = p.posY[j] - posY;
        const float dz = p.posZ[j] - posZ;
        const float r = std::sqrt(dx * dx + dy * dy + dz * dz) + MIN_FLOAT;

        // check for collision distance for all coordinates (instead ternary op)
        if (r > COLLISION_DISTANCE)
        {
            const float divisor = G * dt * p.weight[j] / (r * r * r);
            newVelX += (dx * divisor);
            newVelY += (dy * divisor);
            newVelZ += (dz * divisor);
        }
    }

    tmpVel.x[idx] = newVelX;
    tmpVel.y[idx] = newVelY;
    tmpVel.z[idx] = newVelZ;
} // end of calculate_gravitation_velocity
//----------------------------------------------------------------------------------------------------------------------

/**
 * CUDA kernel to calculate collision velocity
 * @param p      - particles
 * @param tmpVel - temp array for velocities
 * @param N      - Number of particles
 * @param dt     - Size of the time step
 */
__global__ void calculateCollisionVelocity(Particles p, Velocities tmpVel, const unsigned N, float dt)
{
    // calculate thread position in N particles
    const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

    // if out of bounds, end thread work
    if (idx >= N)
        return;

    // obtain thread's particle params ({pos;vel}x,y,z,weight)
    const float posX = p.posX[idx];
    const float posY = p.posY[idx];
    const float posZ = p.posZ[idx];
    const float velX = p.velX[idx];
    const float velY = p.velY[idx];
    const float velZ = p.velZ[idx];
    const float weight = p.weight[idx];

    // place where the new params will be stored
    for (unsigned j = 0u; j < N; ++j)
    {
        const float dx = p.posX[j] - posX;
        const float dy = p.posY[j] - posY;
        const float dz = p.posZ[j] - posZ;
        const float r = std::sqrt(dx * dx + dy * dy + dz * dz);

        if (r > 0.f && r < COLLISION_DISTANCE)
        {
            const float otherWeight = p.weight[j];
            const float divisor = 1 / (weight + otherWeight);
            tmpVel.x[idx] += ((velX * (weight - otherWeight) + 2.f * otherWeight * p.velX[j]) * divisor) - velX;
            tmpVel.y[idx] += ((velY * (weight - otherWeight) + 2.f * otherWeight * p.velY[j]) * divisor) - velY;
            tmpVel.z[idx] += ((velZ * (weight - otherWeight) + 2.f * otherWeight * p.velZ[j]) * divisor) - velZ;
        }
    }

} // end of calculate_collision_velocity
//----------------------------------------------------------------------------------------------------------------------

/**
 * CUDA kernel to update particles
 * @param p      - particles
 * @param tmpVel - temp array for velocities
 * @param N      - Number of particles
 * @param dt     - Size of the time step
 */
__global__ void updateParticles(Particles p, Velocities tmpVel, const unsigned N, float dt)
{
    // calculate thread position in N particles
    const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

    // if out of bounds, end thread work
    if (idx >= N)
        return;

    p.velX[idx] += tmpVel.x[idx];
    p.velY[idx] += tmpVel.y[idx];
    p.velZ[idx] += tmpVel.z[idx];

    p.posX[idx] += p.velX[idx] * dt;
    p.posY[idx] += p.velY[idx] * dt;
    p.posZ[idx] += p.velZ[idx] * dt;
} // end of update_particle
//----------------------------------------------------------------------------------------------------------------------

/**
 * CUDA kernel to calculate particles center of mass
 * @param p    - particles
 * @param com  - pointer to a center of mass
 * @param lock - pointer to a user-implemented lock
 * @param N    - Number of particles
 */
__global__ void centerOfMass(Particles p, float4 *com, int *lock, const unsigned N)
{

} // end of centerOfMass
//----------------------------------------------------------------------------------------------------------------------

/**
 * CPU implementation of the Center of Mass calculation
 * @param particles - All particles in the system
 * @param N         - Number of particles
 */
__host__ float4 centerOfMassRef(MemDesc &memDesc)
{
    float4 com{};

    for (std::size_t i{}; i < memDesc.getDataSize(); i++)
    {
        const float3 pos = {memDesc.getPosX(i), memDesc.getPosY(i), memDesc.getPosZ(i)};
        const float w = memDesc.getWeight(i);

        // Calculate the vector on the line connecting current body and most recent position of center-of-mass
        // Calculate weight ratio only if at least one particle isn't massless
        const float4 d = {pos.x - com.x,
                          pos.y - com.y,
                          pos.z - com.z,
                          ((memDesc.getWeight(i) + com.w) > 0.0f)
                              ? (memDesc.getWeight(i) / (memDesc.getWeight(i) + com.w))
                              : 0.0f};

        // Update position and weight of the center-of-mass according to the weight ration and vector
        com.x += d.x * d.w;
        com.y += d.y * d.w;
        com.z += d.z * d.w;
        com.w += w;
    }

    return com;
} // enf of centerOfMassRef
//----------------------------------------------------------------------------------------------------------------------
