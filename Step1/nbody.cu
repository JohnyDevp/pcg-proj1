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
__global__ void calculateVelocity__(Particles pIn, Particles pOut, const unsigned N, float dt)
{
  // calculate thread position in N particles
  const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
  // if out of bounds, end thread work
  if (idx >= N)
    return;

  // obtain thread's particle params ({pos;vel}x,y,z,weight)
  const float posX = pIn.posX[idx];
  const float posY = pIn.posY[idx];
  const float posZ = pIn.posZ[idx];
  float velX = pIn.velX[idx];
  float velY = pIn.velY[idx];
  float velZ = pIn.velZ[idx];
  const float weight = pIn.weight[idx];

  // place where the new params will be stored
  float newVelX{};
  float newVelY{};
  float newVelZ{};

  const float G_weight = G * weight;
  const float multiplicator = dt / weight;
  // loop through all other particles and compute the differences
  for (unsigned j = 0u; j < N; ++j)
  {
    const float dx = pIn.posX[j] - posX;
    const float dy = pIn.posY[j] - posY;
    const float dz = pIn.posZ[j] - posZ;
    const float r = std::sqrt(dx * dx + dy * dy + dz * dz) + MIN_FLOAT;

    // check for collision distance for all coordinates (instead ternary op)
    if (r > COLLISION_DISTANCE)
    {
      const float divisor = G_weight * pIn.weight[j] * multiplicator / (r * r * r);
      newVelX += (dx * divisor);
      newVelY += (dy * divisor);
      newVelZ += (dz * divisor);
    }
    else if (r > 0.f && r < COLLISION_DISTANCE)
    {
      const float otherWeight = pIn.weight[j];
      const float divisor = 1 / (weight + otherWeight);
      newVelX += ((velX * (weight - otherWeight) + 2.f * otherWeight * pIn.velX[j]) * divisor) - velX;
      newVelY += ((velY * (weight - otherWeight) + 2.f * otherWeight * pIn.velY[j]) * divisor) - velY;
      newVelZ += ((velZ * (weight - otherWeight) + 2.f * otherWeight * pIn.velZ[j]) * divisor) - velZ;
    }
  }

  //= ============================
  velX += newVelX;
  velY += newVelY;
  velZ += newVelZ;

  pOut.posX[idx] = posX + velX * dt;
  pOut.posY[idx] = posY + velY * dt;
  pOut.posZ[idx] = posZ + velZ * dt;

  pOut.velX[idx] = velX;
  pOut.velY[idx] = velY;
  pOut.velZ[idx] = velZ;
}

/**
 * CUDA kernel to calculate new particles velocity and position
 * @param pIn  - particles in
 * @param pOut - particles out
 * @param N    - Number of particles
 * @param dt   - Size of the time step
 */
__global__ void calculateVelocity(Particles pIn, Particles pOut, const unsigned N, float dt)
{
  /***************************************************** DONE *********************************************************/
  /********************************************************************************************************************/
  /*          TODO: CUDA kernel to calculate new particles velocity and position, collapse previous kernels           */
  /********************************************************************************************************************/
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
  float newVelX_gravitation{};
  float newVelY_gravitation{};
  float newVelZ_gravitation{};

  float newVelX_collision{};
  float newVelY_collision{};
  float newVelZ_collision{};
  // loop through all other particles and compute the differences
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
      const float denominator = G * weight * otherWeight / (r * r * r);
      newVelX_gravitation += dx * denominator;
      newVelY_gravitation += dy * denominator;
      newVelZ_gravitation += dz * denominator;
    }
    else if (r < COLLISION_DISTANCE && r > 0.f)
    {
      const float otherWeight = pIn.weight[j];
      const float divisor = 1 / (weight + otherWeight);
      newVelX_collision += ((velX * (weight - otherWeight) + 2.f * otherWeight * pIn.velX[j]) * divisor) - velX;
      newVelY_collision += ((velY * (weight - otherWeight) + 2.f * otherWeight * pIn.velY[j]) * divisor) - velY;
      newVelZ_collision += ((velZ * (weight - otherWeight) + 2.f * otherWeight * pIn.velZ[j]) * divisor) - velZ;
    }
  }
  const float dt_over_weight = dt / weight;
  float newVel_X_Out = newVelX_gravitation * dt_over_weight + newVelX_collision;
  float newVel_Y_Out = newVelY_gravitation * dt_over_weight + newVelY_collision;
  float newVel_Z_Out = newVelZ_gravitation * dt_over_weight + newVelZ_collision;

  //============ update computation

  float velx_final = velX /*original*/ + newVel_X_Out /*new compute vel*/;
  float vely_final = velY /*original*/ + newVel_Y_Out /*new compute vel*/;
  float velz_final = velZ /*original*/ + newVel_Z_Out /*new compute vel*/;

  float posx_final = posX /*original*/ + velx_final * dt;
  float posy_final = posY /*original*/ + vely_final * dt;
  float posz_final = posZ /*original*/ + velz_final * dt;

  // ============ update save
  pOut.posX[idx] = posx_final;
  pOut.posY[idx] = posy_final;
  pOut.posZ[idx] = posz_final;
  pOut.velX[idx] = velx_final;
  pOut.velY[idx] = vely_final;
  pOut.velZ[idx] = velz_final;

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
