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
 * CUDA kernel to calculate new particles velocity and position
 * @param pIn  - particles in
 * @param pOut - particles out
 * @param N    - Number of particles
 * @param dt   - Size of the time step
 */
__global__ void calculateVelocity(Particles pIn, Particles pOut, const unsigned N, float dt)
{
  // calculate thread position in N particles
  unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

  // if out of bounds, end thread work
  if (idx >= N)
    return;

  // params of the current thread's particle
  const float posX = pIn.posX[idx];
  const float posY = pIn.posY[idx];
  const float posZ = pIn.posZ[idx];
  const float velX = pIn.velX[idx];
  const float velY = pIn.velY[idx];
  const float velZ = pIn.velZ[idx];
  const float weight = pIn.weight[idx];

  // place where the new params will be stored
  float newVelX{};
  float newVelY{};
  float newVelZ{};

  // loop through all other particles and compute the differences
  const float G_dt = G * dt;
  for (unsigned j = 0u; j < N; ++j)
  {
    const float dx = pIn.posX[j] - posX;
    const float dy = pIn.posY[j] - posY;
    const float dz = pIn.posZ[j] - posZ;
    const float r = std::sqrt(dx * dx + dy * dy + dz * dz);

    // check for collision distance for all coordinates (instead ternary op)
    if (r > COLLISION_DISTANCE)
    {
      const float otherWeight = pIn.weight[j];
      const float denominator = G_dt * otherWeight / (r * r * r);
      newVelX += dx * denominator;
      newVelY += dy * denominator;
      newVelZ += dz * denominator;
    }
    else if (r < COLLISION_DISTANCE && r > 0.f)
    {
      const float otherWeight = pIn.weight[j];
      const float divisor = 1 / (weight + otherWeight);
      newVelX += ((velX * (weight - otherWeight) + 2.f * otherWeight * pIn.velX[j]) * divisor) - velX;
      newVelY += ((velY * (weight - otherWeight) + 2.f * otherWeight * pIn.velY[j]) * divisor) - velY;
      newVelZ += ((velZ * (weight - otherWeight) + 2.f * otherWeight * pIn.velZ[j]) * divisor) - velZ;
    }
  }

  //============ update computation

  pOut.velX[idx] = velX + newVelX;
  pOut.velY[idx] = velY + newVelY;
  pOut.velZ[idx] = velZ + newVelZ;
  pOut.posX[idx] = posX + (velX + newVelX) * dt;
  pOut.posY[idx] = posY + (velY + newVelY) * dt;
  pOut.posZ[idx] = posZ + (velZ + newVelZ) * dt;

} // end of calculate_gravitation_velocity
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
