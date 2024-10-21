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

__global__ void calculateVelocity(Particles pIn, Particles pOut, const unsigned N, float dt)
{
    extern __shared__ float sharedMem[];

    // calculate thread position in N particles
    const unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

    // Shared memory partitioning
    float* const sharedPosX =     sharedMem;
    float* const sharedPosY =     &sharedPosX[blockDim.x + 1];
    float* const sharedPosZ =     &sharedPosY[blockDim.x + 1];
    float* const sharedVelX =     &sharedPosZ[blockDim.x + 1];
    float* const sharedVelY =     &sharedVelX[blockDim.x + 1];
    float* const sharedVelZ =     &sharedVelY[blockDim.x + 1];
    float* const sharedWeight =   &sharedVelZ[blockDim.x + 1];
    
    // preloaded pointers
    float* const pPosX = pIn.posX;
    float* const pPosY = pIn.posY;
    float* const pPosZ = pIn.posZ;
    float* const pVelX = pIn.velX;
    float* const pVelY = pIn.velY;
    float* const pVelZ = pIn.velZ;
    float* const pWeight = pIn.weight;


    // params of the current thread's particle (if the thread has one)
    const float posX = (idx < N) ? pPosX[idx] : 0.0f;
    const float posY =  (idx < N) ? pPosY[idx] : 0.0f;
    const float posZ =  (idx < N) ? pPosZ[idx] : 0.0f;
    const float velX =  (idx < N) ? pVelX[idx] : 0.0f;
    const float velY =  (idx < N) ? pVelY[idx] : 0.0f;
    const float velZ =  (idx < N) ? pVelZ[idx] : 0.0f;
    const float weight =(idx < N) ? pWeight[idx] : 0.0f;
 

    // place where the new params will be stored
    float newVelX_gravitation{};
    float newVelY_gravitation{};
    float newVelZ_gravitation{};

    float newVelX_collision{};
    float newVelY_collision{};
    float newVelZ_collision{};

    // loop through all other particles and compute the differences
    for (unsigned j = 0u; j < ceilf((float)N / blockDim.x); ++j)
    {
        // load particles to the shared mem
        unsigned loadIdx = j * blockDim.x + threadIdx.x;

        // check whether exists the particle that should be load for current thread
        if (loadIdx < N)
        {
            sharedPosX[threadIdx.x] =   pPosX[loadIdx];
            sharedPosY[threadIdx.x] =   pPosY[loadIdx];
            sharedPosZ[threadIdx.x] =   pPosZ[loadIdx];
            sharedVelX[threadIdx.x] =   pVelX[loadIdx];
            sharedVelY[threadIdx.x] =   pVelY[loadIdx];
            sharedVelZ[threadIdx.x] =   pVelZ[loadIdx];
            sharedWeight[threadIdx.x] = pWeight[loadIdx];
        }

        // sync threads in one block after loading of data
        __syncthreads();

        const int deck = ((j + 1) * blockDim.x > N) ? (N - (j * blockDim.x)) : blockDim.x; //ok
        if (idx > N) // for the threads that serve only for loading
        #   pragma unroll
            for (auto i = 0u; i < deck; i++)
            {

                const float otherPosX = sharedPosX[i];
                const float otherPosY = sharedPosY[i];
                const float otherPosZ = sharedPosZ[i];
                const float otherVelX = sharedVelX[i];
                const float otherVelY = sharedVelY[i];
                const float otherVelZ = sharedVelZ[i];
                const float otherWeight = sharedWeight[i];

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
        
        // sync threads in one block after computation
        __syncthreads();
    }

    if (idx < N)
    {
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
    }

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
