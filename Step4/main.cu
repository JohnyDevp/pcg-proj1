/**
 * @file      main.cu
 *
 * @author    Jan Holáň \n
 *            Faculty of Information Technology \n
 *            Brno University of Technology \n
 *            xholan11@fit.vutbr.cz
 *
 * @brief     PCG Assignment 1
 *
 * @version   2024
 *
 * @date      31 October   2024, 09:00 \n
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
#define CUDA_CALL(call)                                                                                 \
  do                                                                                                    \
  {                                                                                                     \
    const cudaError_t _error = (call);                                                                  \
    if (_error != cudaSuccess)                                                                          \
    {                                                                                                   \
      std::fprintf(stderr, "CUDA error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString(_error)); \
      std::exit(EXIT_FAILURE);                                                                          \
    }                                                                                                   \
  } while (0)

/**
 * Main routine
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
  const unsigned N = static_cast<unsigned>(std::stoul(argv[1]));
  // Length of time step
  const float dt = std::stof(argv[2]);
  // Number of steps
  const unsigned steps = static_cast<unsigned>(std::stoul(argv[3]));
  // Number of thread blocks
  const unsigned simBlockDim = static_cast<unsigned>(std::stoul(argv[4]));
  // Write frequency
  const unsigned writeFreq = static_cast<unsigned>(std::stoul(argv[5]));
  // number of reduction threads
  const unsigned redTotalThreadCount = static_cast<unsigned>(std::stoul(argv[6]));
  // Number of reduction threads/blocks
  const unsigned redBlockDim = static_cast<unsigned>(std::stoul(argv[7]));

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
  float4 *hCenterOfMass{};

  /********************************************************************************************************************/
  /*                                    CPU side memory allocation (pinned)                                           */
  /********************************************************************************************************************/

  size_t size = sizeof(float) * N;

  CUDA_CALL(cudaHostAlloc(&hParticles.posX, size, cudaHostAllocDefault));
  CUDA_CALL(cudaHostAlloc(&hParticles.posY, size, cudaHostAllocDefault));
  CUDA_CALL(cudaHostAlloc(&hParticles.posZ, size, cudaHostAllocDefault));
  CUDA_CALL(cudaHostAlloc(&hParticles.velX, size, cudaHostAllocDefault));
  CUDA_CALL(cudaHostAlloc(&hParticles.velY, size, cudaHostAllocDefault));
  CUDA_CALL(cudaHostAlloc(&hParticles.velZ, size, cudaHostAllocDefault));
  CUDA_CALL(cudaHostAlloc(&hParticles.weight, size, cudaHostAllocDefault));

  CUDA_CALL(cudaHostAlloc(&hCenterOfMass, sizeof(float4), cudaHostAllocDefault));

  /********************************************************************************************************************/
  /*                                    Fill memory descriptor layout                                                 */
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
  catch (const std::exception &e)
  {
    std::fprintf(stderr, "Error: %s\n", e.what());
    return EXIT_FAILURE;
  }

  Particles dParticles[2]{};
  float4 *dCenterOfMass{};
  int *dLock{};

  /********************************************************************************************************************/
  /*                                           GPU side memory allocation                                             */
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

  CUDA_CALL(cudaMalloc<float4>(&dCenterOfMass, sizeof(float4)));
  CUDA_CALL(cudaMalloc<int>(&dLock, sizeof(int)));

  /********************************************************************************************************************/
  /*                                           Memory transfer CPU -> GPU                                             */
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

  /********************************************************************************************************************/
  /*                                           Clear GPU center of mass                                               */
  /********************************************************************************************************************/
  CUDA_CALL(cudaMemset(dCenterOfMass, 0, sizeof(float4)));
  CUDA_CALL(cudaMemset(dLock, 0, sizeof(int)));

  /********************************************************************************************************************/
  /*                           TODO: Declare and create necessary CUDA streams and events                             */
  /********************************************************************************************************************/
  cudaStream_t computeVelocityStream, transferStream, computeMassStream;
  CUDA_CALL(cudaStreamCreate(&computeVelocityStream));
  CUDA_CALL(cudaStreamCreate(&transferStream));
  CUDA_CALL(cudaStreamCreate(&computeMassStream));

  // Create CUDA events for synchronization
  cudaEvent_t computeVelocityEvent, transferEvent, computeMassEvent;
  CUDA_CALL(cudaEventCreate(&computeVelocityEvent));
  CUDA_CALL(cudaEventCreate(&computeMassEvent));
  CUDA_CALL(cudaEventCreate(&transferEvent));

  // Get CUDA device warp size
  int device;
  int warpSize;

  CUDA_CALL(cudaGetDevice(&device));
  CUDA_CALL(cudaDeviceGetAttribute(&warpSize, cudaDevAttrWarpSize, device));

  /********************************************************************************************************************/
  /*                                        Set dynamic shared memory computation                                     */
  /********************************************************************************************************************/
  const std::size_t sharedMemSize = simBlockDim * (7 * sizeof(float)) + 1;
  const std::size_t redSharedMemSize = redBlockDim * sizeof(float4); // you can use warpSize variable

  // Lambda for checking if we should write current step to the file
  auto shouldWrite = [writeFreq](unsigned s) -> bool
  {
    return writeFreq > 0u && (s % writeFreq == 0u);
  };

  // Lamda for getting record number
  auto getRecordNum = [writeFreq](unsigned s) -> unsigned
  {
    return s / writeFreq;
  };

  // Start measurement
  const auto start = std::chrono::steady_clock::now();

  /********************************************************************************************************************/
  /*            TODO: Edit the loop to work asynchronously and overlap computation with data transfers.               */
  /*                  Use shouldWrite lambda to determine if data should be outputted to file.                        */
  /*                           if (shouldWrite(s, writeFreq)) { ... }                                                 */
  /*                        Use getRecordNum lambda to get the record number.                                         */
  /********************************************************************************************************************/

  for (unsigned s = 0u; s < steps; ++s)
  {
    const unsigned srcIdx = s % 2;       // source particles index
    const unsigned dstIdx = (s + 1) % 2; // destination particles index

    // kernel for position update in compute stream
    calculateVelocity<<<simGridDim, simBlockDim, sharedMemSize, computeVelocityStream>>>(dParticles[srcIdx], dParticles[dstIdx], N, dt);
    CUDA_CALL(cudaEventRecord(computeVelocityEvent, computeVelocityStream));

    if (shouldWrite(s))
    {
      auto recordNum = getRecordNum(s);
      // prevent memory goes before previous kernel finishes
      CUDA_CALL(cudaStreamWaitEvent(transferStream, computeVelocityEvent));

      // transfer particle data back to the CPU asynchronously
      CUDA_CALL(cudaMemcpyAsync(hParticles.posX, dParticles[srcIdx].posX, N * sizeof(float), cudaMemcpyDeviceToHost, transferStream));
      CUDA_CALL(cudaMemcpyAsync(hParticles.posY, dParticles[srcIdx].posY, N * sizeof(float), cudaMemcpyDeviceToHost, transferStream));
      CUDA_CALL(cudaMemcpyAsync(hParticles.posZ, dParticles[srcIdx].posZ, N * sizeof(float), cudaMemcpyDeviceToHost, transferStream));
      CUDA_CALL(cudaMemcpyAsync(hParticles.velX, dParticles[srcIdx].velX, N * sizeof(float), cudaMemcpyDeviceToHost, transferStream));
      CUDA_CALL(cudaMemcpyAsync(hParticles.velY, dParticles[srcIdx].velY, N * sizeof(float), cudaMemcpyDeviceToHost, transferStream));
      CUDA_CALL(cudaMemcpyAsync(hParticles.velZ, dParticles[srcIdx].velZ, N * sizeof(float), cudaMemcpyDeviceToHost, transferStream));
      CUDA_CALL(cudaMemcpyAsync(hParticles.weight, dParticles[srcIdx].weight, N * sizeof(float), cudaMemcpyDeviceToHost, transferStream));

      CUDA_CALL(cudaEventRecord(transferEvent, transferStream));
      CUDA_CALL(cudaEventSynchronize(transferEvent));
      h5Helper.writeParticleData(recordNum);

      // zero out center of mass and compute centerOfMass
      CUDA_CALL(cudaMemsetAsync(dCenterOfMass, 0, sizeof(float4), computeMassStream));
      centerOfMass<<<redGridDim, redBlockDim, redSharedMemSize, computeMassStream>>>(dParticles[srcIdx], dCenterOfMass, dLock, N);

      // record an event when center of mass calculation is done
      CUDA_CALL(cudaEventRecord(computeMassEvent, computeMassStream));

      // wait for the center of mass calculation to complete before transfer
      CUDA_CALL(cudaStreamWaitEvent(transferStream, computeMassEvent));

      // transfer center of mass to CPU and then write it
      CUDA_CALL(cudaMemcpyAsync(hCenterOfMass, dCenterOfMass, sizeof(float4), cudaMemcpyDeviceToHost, transferStream));

      // now all memory events are recorded, wait for all transfers to finish and write the data
      CUDA_CALL(cudaEventRecord(transferEvent, transferStream));
      CUDA_CALL(cudaEventSynchronize(transferEvent));
      h5Helper.writeCom(*hCenterOfMass, recordNum);
    }
  }

const unsigned resIdx = steps % 2; // result particles index

  /********************************************************************************************************************/
  /*                          TODO: Invocation of center of mass kernel, do not forget to add                         */
  /*                              additional synchronization and set appropriate stream                               */
  /********************************************************************************************************************/
  CUDA_CALL(cudaStreamWaitEvent(computeMassStream, computeVelocityEvent));
  CUDA_CALL(cudaMemsetAsync(dCenterOfMass, 0, sizeof(float4), computeMassStream));
  centerOfMass<<<redGridDim, redBlockDim, redSharedMemSize, computeMassStream>>>(dParticles[resIdx], dCenterOfMass, dLock, N);

  // ================================================================================================================== //

  // Wait for all CUDA kernels to finish
  CUDA_CALL(cudaDeviceSynchronize());

  // End measurement
  const auto end = std::chrono::steady_clock::now();

  // Approximate simulation wall time
  const float elapsedTime = std::chrono::duration<float>(end - start).count();
  std::printf("Time: %f s\n", elapsedTime);

  /********************************************************************************************************************/
  /*                                           Memory transfer GPU -> CPU                                             */
  /********************************************************************************************************************/
  CUDA_CALL(cudaMemcpy(hParticles.posX, dParticles[resIdx].posX, N * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(hParticles.posY, dParticles[resIdx].posY, N * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(hParticles.posZ, dParticles[resIdx].posZ, N * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(hParticles.velX, dParticles[resIdx].velX, N * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(hParticles.velY, dParticles[resIdx].velY, N * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(hParticles.velZ, dParticles[resIdx].velZ, N * sizeof(float), cudaMemcpyDeviceToHost));
  CUDA_CALL(cudaMemcpy(hParticles.weight, dParticles[resIdx].weight, N * sizeof(float), cudaMemcpyDeviceToHost));

  CUDA_CALL(cudaMemcpy(hCenterOfMass, dCenterOfMass, sizeof(float4), cudaMemcpyDeviceToHost));

  // Compute reference center of mass on CPU
  const float4 refCenterOfMass = centerOfMassRef(md);

  std::printf("Reference center of mass: %f, %f, %f, %f\n",
              refCenterOfMass.x,
              refCenterOfMass.y,
              refCenterOfMass.z,
              refCenterOfMass.w);

  std::printf("Center of mass on GPU: %f, %f, %f, %f\n",
              hCenterOfMass->x,
              hCenterOfMass->y,
              hCenterOfMass->z,
              hCenterOfMass->w);

  // Writing final values to the file
  h5Helper.writeComFinal(*hCenterOfMass);
  h5Helper.writeParticleDataFinal();

  /********************************************************************************************************************/
  /*                                  TODO: CUDA streams and events destruction                                       */
  /********************************************************************************************************************/
  CUDA_CALL(cudaStreamDestroy(computeVelocityStream));
  CUDA_CALL(cudaStreamDestroy(transferStream));
  CUDA_CALL(cudaStreamDestroy(computeMassStream));

  CUDA_CALL(cudaEventDestroy(computeVelocityEvent));
  CUDA_CALL(cudaEventDestroy(computeMassEvent));
  CUDA_CALL(cudaEventDestroy(transferEvent));

  /***************************************************** DONE *********************************************************/
  /********************************************************************************************************************/
  /*                                           GPU side memory deallocation                                           */
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
  CUDA_CALL(cudaFree(dCenterOfMass));
  CUDA_CALL(cudaFree(dLock));

  /***************************************************** DONE *********************************************************/
  /********************************************************************************************************************/
  /*                                           CPU side memory deallocation                                           */
  /********************************************************************************************************************/
  CUDA_CALL(cudaFreeHost(hParticles.posX));
  CUDA_CALL(cudaFreeHost(hParticles.posY));
  CUDA_CALL(cudaFreeHost(hParticles.posZ));
  CUDA_CALL(cudaFreeHost(hParticles.velX));
  CUDA_CALL(cudaFreeHost(hParticles.velY));
  CUDA_CALL(cudaFreeHost(hParticles.velZ));
  CUDA_CALL(cudaFreeHost(hParticles.weight));
  CUDA_CALL(cudaFreeHost(hCenterOfMass));

} // end of main
//----------------------------------------------------------------------------------------------------------------------
