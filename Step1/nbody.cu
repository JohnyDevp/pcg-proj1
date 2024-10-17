/**
 * @file      nbody.cu
 *
 * @author    Name Surname \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            xlogin00@fit.vutbr.cz
 *
 * @brief     PCG Assignment 1
 *
 * @version   2024
 *
 * @date      04 October   2023, 09:00 (created) \n
 */

#include <device_launch_parameters.h>

#include "nbody.cuh"

/* Constants */
constexpr float G                  = 6.67384e-11f;
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
  /***************************************************** DONE *********************************************************/
  /********************************************************************************************************************/
  /*          TODO: CUDA kernel to calculate new particles velocity and position, collapse previous kernels           */
  /********************************************************************************************************************/
    // preloaded pointers
    float* const pPosX = pIn.posX;
    float* const pPosY = pIn.posY;
    float* const pPosZ = pIn.posZ;
    float* const pVelX = pIn.velX;
    float* const pVelY = pIn.velY;
    float* const pVelZ = pIn.velZ;
    float* const pWeight = pIn.weight;

    // calculate thread position in N particles
    unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

    // if out of bounds, end thread work
    if (idx >= N) return;

    // params of the current thread's particle
    const float posX = pPosX[idx];
    const float posY = pPosY[idx];
    const float posZ = pPosZ[idx];
    const float velX = pVelX[idx];
    const float velY = pVelY[idx];
    const float velZ = pVelZ[idx];
    const float weight = pWeight[idx];

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
        const float otherPosX = pPosX[j];
        const float otherPosY = pPosY[j];
        const float otherPosZ = pPosZ[j];
        const float otherVelX = pVelX[j];
        const float otherVelY = pVelY[j];
        const float otherVelZ = pVelZ[j];
        const float otherWeight = pWeight[j];

        const float dx = otherPosX - posX;
        const float dy = otherPosY - posY;
        const float dz = otherPosZ - posZ;

        const float r2 = dx * dx + dy * dy + dz * dz;
        const float r = std::sqrt(r2);
        const float f = G * weight * otherWeight / r2 + MIN_FLOAT;

        // check for collision distance for all coordinates (instead ternary op)
        if ((r + MIN_FLOAT) > COLLISION_DISTANCE) {
            const float denominator = f / r;
            newVelX_gravitation += dx * denominator;
            newVelY_gravitation += dy * denominator;
            newVelZ_gravitation += dz * denominator;
        } 
        else if (r < COLLISION_DISTANCE && r > 0.f)
        {
            newVelX_collision += ((weight * velX - otherWeight * velX + 2.f * otherWeight * otherVelX) / (weight + otherWeight)) - velX;
            newVelY_collision += ((weight * velY - otherWeight * velY + 2.f * otherWeight * otherVelY) / (weight + otherWeight)) - velY;
            newVelZ_collision += ((weight * velZ - otherWeight * velZ + 2.f * otherWeight * otherVelZ) / (weight + otherWeight)) - velZ;
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

}// end of calculate_gravitation_velocity
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
