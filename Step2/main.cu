/**
 * @file      main.cu
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

#include <cmath>
#include <cstdio>
#include <chrono>
#include <string>

#include "nbody.cuh"
#include "h5Helper.h"

/**
 * @brief CUDA error checking macro
 * @param call CUDA API call
 */
#define CUDA_CALL(call) \
  do { \
    const cudaError_t _error = (call); \
    if (_error != cudaSuccess) \
    { \
      std::fprintf(stderr, "CUDA error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(_error)); \
      std::exit(EXIT_FAILURE); \
    } \
  } while(0)

/**
 * Main rotine
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv)
{
  if (argc != 10)
  {
    std::printf("Usage: nbody <N> <dt> <steps> <threads/block> <write intesity> <reduction threads> <reduction threads/block> <input> <output>\n");
    std::exit(1);
  }

  // Number of particles
  const unsigned N                   = static_cast<unsigned>(std::stoul(argv[1]));
  // Length of time step
  const float    dt                  = std::stof(argv[2]);
  // Number of steps
  const unsigned steps               = static_cast<unsigned>(std::stoul(argv[3]));
  // Number of thread blocks
  const unsigned simBlockDim         = static_cast<unsigned>(std::stoul(argv[4]));
  // Write frequency
  const unsigned writeFreq           = static_cast<unsigned>(std::stoul(argv[5]));
  // number of reduction threads
  const unsigned redTotalThreadCount = static_cast<unsigned>(std::stoul(argv[6]));
  // Number of reduction threads/blocks
  const unsigned redBlockDim         = static_cast<unsigned>(std::stoul(argv[7]));

  // Size of the simulation CUDA grid - number of blocks
  const unsigned simGridDim = (N + simBlockDim - 1) / simBlockDim;
  // Size of the reduction CUDA grid - number of blocks
  const unsigned redGridDim = (redTotalThreadCount + redBlockDim - 1) / redBlockDim;

  // Log benchmark setup
  std::printf("       NBODY GPU simulation\n"
              "N:                       %u\n"
              "dt:                      %f\n"
              "steps:                   %u\n"
              "threads/block:           %u\n"
              "blocks/grid:             %u\n"
              "reduction threads/block: %u\n"
              "reduction blocks/grid:   %u\n",
              N, dt, steps, simBlockDim, simGridDim, redBlockDim, redGridDim);

  const std::size_t recordsCount = (writeFreq > 0) ? (steps + writeFreq - 1) / writeFreq : 0;

  Particles hParticles{};

  /***************************************************** DONE *********************************************************/
  /********************************************************************************************************************/
  /*                              TODO: CPU side memory allocation (pinned)                                           */
  /********************************************************************************************************************/

  size_t size = sizeof(float) * N;

  CUDA_CALL(cudaHostAlloc(&hParticles.posX, size, cudaHostAllocDefault));
  CUDA_CALL(cudaHostAlloc(&hParticles.posY, size, cudaHostAllocDefault));
  CUDA_CALL(cudaHostAlloc(&hParticles.posZ, size, cudaHostAllocDefault));
  CUDA_CALL(cudaHostAlloc(&hParticles.velX, size, cudaHostAllocDefault));
  CUDA_CALL(cudaHostAlloc(&hParticles.velY, size, cudaHostAllocDefault));
  CUDA_CALL(cudaHostAlloc(&hParticles.velZ, size, cudaHostAllocDefault));
  CUDA_CALL(cudaHostAlloc(&hParticles.weight, size, cudaHostAllocDefault));

  /***************************************************** DONE *********************************************************/
  /********************************************************************************************************************/
  /*                              TODO: Fill memory descriptor layout                                                 */
  /********************************************************************************************************************/
  /*
   * Caution! Create only after CPU side allocation
   * parameters:
   *                            Stride of two            Offset of the first
   *       Data pointer       consecutive elements        element in FLOATS,
   *                          in FLOATS, not bytes            not bytes
  */
  MemDesc md(hParticles.posX, 1, 0,
      hParticles.posY, 1, 0,
      hParticles.posZ, 1, 0,
      hParticles.velX, 1, 0,
      hParticles.velY, 1, 0,
      hParticles.velZ, 1, 0,
      hParticles.weight, 1, 0,
      N,
      recordsCount);

  // Initialisation of helper class and loading of input data
  H5Helper h5Helper(argv[8], argv[9], md);

  try
  {
    h5Helper.init();
    h5Helper.readParticleData();
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "Error: %s\n", e.what());
    return EXIT_FAILURE;
  }

  Particles dParticles[2]{};

  /***************************************************** DONE *********************************************************/
  /********************************************************************************************************************/
  /*                                     TODO: GPU side memory allocation                                             */
  /********************************************************************************************************************/

  for (auto i = 0u; i < 2; i++)
  {
      CUDA_CALL(cudaMalloc<float>(&dParticles[i].posX, size));
      CUDA_CALL(cudaMalloc<float>(&dParticles[i].posY, size));
      CUDA_CALL(cudaMalloc<float>(&dParticles[i].posZ, size));
      CUDA_CALL(cudaMalloc<float>(&dParticles[i].velX, size));
      CUDA_CALL(cudaMalloc<float>(&dParticles[i].velY, size));
      CUDA_CALL(cudaMalloc<float>(&dParticles[i].velZ, size));
      CUDA_CALL(cudaMalloc<float>(&dParticles[i].weight, size));
  }
  
  /***************************************************** DONE *********************************************************/
  /********************************************************************************************************************/
  /*                                     TODO: Memory transfer CPU -> GPU                                             */
  /********************************************************************************************************************/
  for (auto i = 0u; i < 2; i++)
  {
      CUDA_CALL(cudaMemcpy(dParticles[i].posX, hParticles.posX, N * sizeof(float), cudaMemcpyHostToDevice));
      CUDA_CALL(cudaMemcpy(dParticles[i].posY, hParticles.posY, N * sizeof(float), cudaMemcpyHostToDevice));
      CUDA_CALL(cudaMemcpy(dParticles[i].posZ, hParticles.posZ, N * sizeof(float), cudaMemcpyHostToDevice));
      CUDA_CALL(cudaMemcpy(dParticles[i].velX, hParticles.velX, N * sizeof(float), cudaMemcpyHostToDevice));
      CUDA_CALL(cudaMemcpy(dParticles[i].velY, hParticles.velY, N * sizeof(float), cudaMemcpyHostToDevice));
      CUDA_CALL(cudaMemcpy(dParticles[i].velZ, hParticles.velZ, N * sizeof(float), cudaMemcpyHostToDevice));
      CUDA_CALL(cudaMemcpy(dParticles[i].weight, hParticles.weight, N * sizeof(float), cudaMemcpyHostToDevice));
  }

  /***************************************************** DONE *********************************************************/
  /********************************************************************************************************************/
  /*                                  TODO: Set dynamic shared memory computation                                     */
  /********************************************************************************************************************/
  const std::size_t sharedMemSize = simBlockDim * (7 * sizeof(float)) + 1;

  // Start measurement
  const auto start = std::chrono::steady_clock::now();

  for (unsigned s = 0u; s < steps; ++s)
  {
    const unsigned srcIdx = s % 2;        // source particles index
    const unsigned dstIdx = (s + 1) % 2;  // destination particles index

    /***************************************************** DONE *******************************************************/
    /******************************************************************************************************************/
    /*                   TODO: GPU kernel invocation with correctly set dynamic memory size                           */
    /******************************************************************************************************************/
    calculateVelocity<<<simGridDim, simBlockDim, sharedMemSize>>>(dParticles[srcIdx], dParticles[dstIdx], N, dt);

  }

  // Wait for all CUDA kernels to finish
  CUDA_CALL(cudaDeviceSynchronize());

  // End measurement
  const auto end = std::chrono::steady_clock::now();

  // Approximate simulation wall time
  const float elapsedTime = std::chrono::duration<float>(end - start).count();
  std::printf("Time: %f s\n", elapsedTime);

  const unsigned resIdx = steps % 2;    // result particles index

  /***************************************************** DONE *********************************************************/
  /********************************************************************************************************************/
  /*                                     TODO: Memory transfer GPU -> CPU                                             */
  /********************************************************************************************************************/
  CUDA_CALL(cudaMemcpy(hParticles.posX, dParticles[resIdx].posX, N * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(hParticles.posY, dParticles[resIdx].posY, N * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(hParticles.posZ, dParticles[resIdx].posZ, N * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(hParticles.velX, dParticles[resIdx].velX, N * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(hParticles.velY, dParticles[resIdx].velY, N * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(hParticles.velZ, dParticles[resIdx].velZ, N * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(hParticles.weight, dParticles[resIdx].weight, N * sizeof(float), cudaMemcpyDeviceToHost));



  // Compute reference center of mass on CPU
  const float4 refCenterOfMass = centerOfMassRef(md);

  std::printf("Reference center of mass: %f, %f, %f, %f\n",
              refCenterOfMass.x,
              refCenterOfMass.y,
              refCenterOfMass.z,
              refCenterOfMass.w);

  std::printf("Center of mass on GPU: %f, %f, %f, %f\n", 0.f, 0.f, 0.f, 0.f);

  // Writing final values to the file
  h5Helper.writeComFinal(refCenterOfMass);
  h5Helper.writeParticleDataFinal();

  /***************************************************** DONE *********************************************************/
 /********************************************************************************************************************/
 /*                                     TODO: GPU side memory deallocation                                           */
 /********************************************************************************************************************/

  for (auto i = 0u; i < 2; i++)
  {
      CUDA_CALL(cudaFree(dParticles[i].posX));
      CUDA_CALL(cudaFree(dParticles[i].posY));
      CUDA_CALL(cudaFree(dParticles[i].posZ));
      CUDA_CALL(cudaFree(dParticles[i].velX));
      CUDA_CALL(cudaFree(dParticles[i].velY));
      CUDA_CALL(cudaFree(dParticles[i].velZ));
      CUDA_CALL(cudaFree(dParticles[i].weight));
  }

  /***************************************************** DONE *********************************************************/
  /********************************************************************************************************************/
  /*                                     TODO: CPU side memory deallocation                                           */
  /********************************************************************************************************************/

  CUDA_CALL(cudaFreeHost(hParticles.posX));
  CUDA_CALL(cudaFreeHost(hParticles.posY));
  CUDA_CALL(cudaFreeHost(hParticles.posZ));
  CUDA_CALL(cudaFreeHost(hParticles.velX));
  CUDA_CALL(cudaFreeHost(hParticles.velY));
  CUDA_CALL(cudaFreeHost(hParticles.velZ));
  CUDA_CALL(cudaFreeHost(hParticles.weight));


}// end of main
//----------------------------------------------------------------------------------------------------------------------
