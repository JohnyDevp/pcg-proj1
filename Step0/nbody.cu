/**
 * @file      nbody.cu
 *
 * @author    Jan Hol·Ú \n
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
constexpr float G                  = 6.67384e-11f;
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
    /***************************************************** DONE *********************************************************/
    /********************************************************************************************************************/
    /*              TODO: CUDA kernel to calculate gravitation velocity, see reference CPU version                      */
    /********************************************************************************************************************/

    // preloaded pointers
    float* const pPosX = p.posX;
    float* const pPosY = p.posY;
    float* const pPosZ = p.posZ;
    float* const pVelX = p.velX;
    float* const pVelY = p.velY;
    float* const pVelZ = p.velZ;
    float* const pWeight = p.weight;

    // calculate thread position in N particles
    unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

    // if out of bounds, end thread work
    if (idx >= N) return;

    // params of the current thread's particle
    const float posX = pPosX[idx];
    const float posY = pPosY[idx];
    const float posZ = pPosZ[idx];
    const float weight = pWeight[idx];

    // place where the new params will be stored
    float newVelX{};
    float newVelY{};
    float newVelZ{};

    // loop through all other particles and compute the differences
    for (unsigned j = 0u; j < N; ++j)
    {
        const float otherPosX = pPosX[j];
        const float otherPosY = pPosY[j];
        const float otherPosZ = pPosZ[j];
        const float otherWeight = pWeight[j];

        const float dx = otherPosX - posX;
        const float dy = otherPosY - posY;
        const float dz = otherPosZ - posZ;

        const float r2 = dx * dx + dy * dy + dz * dz;
        const float r = std::sqrt(r2) + MIN_FLOAT;

        const float f = G * weight * otherWeight / r2 + MIN_FLOAT;

        // check for collision distance for all coordinates (instead ternary op)
        if (r > COLLISION_DISTANCE) {
            const float denominator = f / r;
            newVelX += dx * denominator;
            newVelY += dy * denominator;
            newVelZ += dz * denominator;
        }
    }

    tmpVel.x[idx] = newVelX * (dt / weight);
    tmpVel.y[idx] = newVelY * (dt / weight);
    tmpVel.z[idx] = newVelZ * (dt / weight);
}// end of calculate_gravitation_velocity
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
    /***************************************************** DONE *********************************************************/
    /********************************************************************************************************************/
    /*              TODO: CUDA kernel to calculate collision velocity, see reference CPU version                        */
    /********************************************************************************************************************/
   
    // calculate thread position in N particles
    unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

    // if out of bounds, end thread work
    if (idx >= N) return;
    
    // obtain thread's particle params ({pos;vel}x,y,z,weight)
    const float posX = p.posX[idx];
    const float posY = p.posY[idx];
    const float posZ = p.posZ[idx];
    const float velX = p.velX[idx];
    const float velY = p.velY[idx];
    const float velZ = p.velZ[idx];
    const float weight = p.weight[idx];

    // place where the new params will be stored
    float newVelX{};
    float newVelY{};
    float newVelZ{};

    for (unsigned j = 0u; j < N; ++j)
    {
        const float otherPosX = p.posX[j];
        const float otherPosY = p.posY[j];
        const float otherPosZ = p.posZ[j];
        const float otherVelX = p.velX[j];
        const float otherVelY = p.velY[j];
        const float otherVelZ = p.velZ[j];
        const float otherWeight = p.weight[j];

        const float dx = otherPosX - posX;
        const float dy = otherPosY - posY;
        const float dz = otherPosZ - posZ;

        const float r2 = dx * dx + dy * dy + dz * dz;
        const float r = std::sqrt(r2);

        if (r < COLLISION_DISTANCE && r > 0.f)
        {
            newVelX += ((weight * velX - otherWeight * velX + 2.f * otherWeight * otherVelX) / (weight + otherWeight)) - velX;
            newVelY += ((weight * velY - otherWeight * velY + 2.f * otherWeight * otherVelY) / (weight + otherWeight)) - velY;
            newVelZ += ((weight * velZ - otherWeight * velZ + 2.f * otherWeight * otherVelZ) / (weight + otherWeight)) - velZ;
        }
    }

    tmpVel.x[idx] += newVelX;
    tmpVel.y[idx] += newVelY;
    tmpVel.z[idx] += newVelZ;
    

}// end of calculate_collision_velocity
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
    /***************************************************** DONE *********************************************************/
    /********************************************************************************************************************/
    /*             TODO: CUDA kernel to update particles velocities and positions, see reference CPU version            */
    /********************************************************************************************************************/

    // calculate thread position in N particles
    unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

    // if out of bounds, end thread work
    if (idx >= N) return;

    // take the pointers as shortcut per thread
    float* const pPosX = p.posX;
    float* const pPosY = p.posY;
    float* const pPosZ = p.posZ;
    float* const pVelX = p.velX;
    float* const pVelY = p.velY;
    float* const pVelZ = p.velZ;
    float* const pWeight = p.weight;

    float* const tmpVelX = tmpVel.x;
    float* const tmpVelY = tmpVel.y;
    float* const tmpVelZ = tmpVel.z;

    // load current values
    float posX = pPosX[idx];
    float posY = pPosY[idx];
    float posZ = pPosZ[idx];

    float velX = pVelX[idx];
    float velY = pVelY[idx];
    float velZ = pVelZ[idx];


    // compute new values
    velX += tmpVelX[idx];
    velY += tmpVelY[idx];
    velZ += tmpVelZ[idx];

    posX += velX * dt;
    posY += velY * dt;
    posZ += velZ * dt;

    // save new values
    pPosX[idx] = posX;
    pPosY[idx] = posY;
    pPosZ[idx] = posZ;

    pVelX[idx] = velX;
    pVelY[idx] = velY;
    pVelZ[idx] = velZ;
    

}// end of update_particle
//----------------------------------------------------------------------------------------------------------------------

/**
 * CUDA kernel to calculate particles center of mass
 * @param p    - particles
 * @param com  - pointer to a center of mass
 * @param lock - pointer to a user-implemented lock
 * @param N    - Number of particles
 */
__global__ void centerOfMass(Particles p, float4* com, int* lock, const unsigned N)
{

}// end of centerOfMass
//----------------------------------------------------------------------------------------------------------------------

/**
 * CPU implementation of the Center of Mass calculation
 * @param particles - All particles in the system
 * @param N         - Number of particles
 */
__host__ float4 centerOfMassRef(MemDesc& memDesc)
{
  float4 com{};

  for (std::size_t i{}; i < memDesc.getDataSize(); i++)
  {
    const float3 pos = {memDesc.getPosX(i), memDesc.getPosY(i), memDesc.getPosZ(i)};
    const float  w   = memDesc.getWeight(i);

    // Calculate the vector on the line connecting current body and most recent position of center-of-mass
    // Calculate weight ratio only if at least one particle isn't massless
    const float4 d = {pos.x - com.x,
                      pos.y - com.y,
                      pos.z - com.z,
                      ((memDesc.getWeight(i) + com.w) > 0.0f)
                        ? ( memDesc.getWeight(i) / (memDesc.getWeight(i) + com.w))
                        : 0.0f};

    // Update position and weight of the center-of-mass according to the weight ration and vector
    com.x += d.x * d.w;
    com.y += d.y * d.w;
    com.z += d.z * d.w;
    com.w += w;
  }

  return com;
}// enf of centerOfMassRef
//----------------------------------------------------------------------------------------------------------------------
