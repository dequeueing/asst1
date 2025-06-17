#include <stdio.h>
#include <thread>
#include <algorithm>

#include "CycleTimer.h"

// Forward declaration of the mandel function
static inline int mandel(float c_re, float c_im, int count);

typedef struct {
    float x0, x1;
    float y0, y1;
    unsigned int width;
    unsigned int height;
    int maxIterations;
    int* output;
    int threadId;
    int numThreads;
} WorkerArgs;


extern void mandelbrotSerial(
    float x0, float y0, float x1, float y1,
    int width, int height,
    int startRow, int numRows,
    int maxIterations,
    int output[]);


//
// workerThreadStart --
//
// Thread entrypoint.
void workerThreadStart(WorkerArgs * const args) {
    // Start timing
    double startTime = CycleTimer::currentSeconds();
    
    // Calculate coordinate deltas (same as serial implementation)
    float dx = (args->x1 - args->x0) / args->width;
    float dy = (args->y1 - args->y0) / args->height;
    
    // Define block size (32x32 pixels)
    const int BLOCK_SIZE = 32;
    
    // Calculate number of blocks in each dimension
    int blocksX = (args->width + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int blocksY = (args->height + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int totalBlocks = blocksX * blocksY;
    
    // Each thread takes every nth block where n is the number of threads
    for (int block = args->threadId; block < totalBlocks; block += args->numThreads) {
        int blockX = block % blocksX;
        int blockY = block / blocksX;
        
        // Calculate start and end coordinates for this block
        int startX = blockX * BLOCK_SIZE;
        int startY = blockY * BLOCK_SIZE;
        int endX = std::min<int>(startX + BLOCK_SIZE, args->width);
        int endY = std::min<int>(startY + BLOCK_SIZE, args->height);
        
        // Process each pixel in this block
        for (int y = startY; y < endY; y++) {
            for (int x = startX; x < endX; x++) {
                // Use same coordinate calculation as serial version
                float x_coord = args->x0 + x * dx;
                float y_coord = args->y0 + y * dy;
                int index = y * args->width + x;
                args->output[index] = mandel(x_coord, y_coord, args->maxIterations);
            }
        }
    }
    
    // End timing
    double endTime = CycleTimer::currentSeconds();
    double threadTime = endTime - startTime;
    
    // Print timing information
    // printf("Thread %d: processed %d blocks in %.3f ms\n",
    //        args->threadId, (totalBlocks - args->threadId + args->numThreads - 1) / args->numThreads,
    //        threadTime * 1000);
}

//
// MandelbrotThread --
//
// Multi-threaded implementation of mandelbrot set image generation.
// Threads of execution are created by spawning std::threads.
void mandelbrotThread(
    int numThreads,
    float x0, float y0, float x1, float y1,
    int width, int height,
    int maxIterations, int output[])
{
    static constexpr int MAX_THREADS = 32;

    if (numThreads > MAX_THREADS)
    {
        fprintf(stderr, "Error: Max allowed threads is %d\n", MAX_THREADS);
        exit(1);
    }

    // Creates thread objects that do not yet represent a thread.
    std::thread workers[MAX_THREADS];
    WorkerArgs args[MAX_THREADS];

    for (int i=0; i<numThreads; i++) {
      
        // TODO FOR CS149 STUDENTS: You may or may not wish to modify
        // the per-thread arguments here.  The code below copies the
        // same arguments for each thread
        args[i].x0 = x0;
        args[i].y0 = y0;
        args[i].x1 = x1;
        args[i].y1 = y1;
        args[i].width = width;
        args[i].height = height;
        args[i].maxIterations = maxIterations;
        args[i].numThreads = numThreads;
        args[i].output = output;
      
        args[i].threadId = i;
    }

    // Spawn the worker threads.  Note that only numThreads-1 std::threads
    // are created and the main application thread is used as a worker
    // as well.
    for (int i=1; i<numThreads; i++) {
        workers[i] = std::thread(workerThreadStart, &args[i]);
    }
    
    workerThreadStart(&args[0]);

    // join worker threads
    for (int i=1; i<numThreads; i++) {
        workers[i].join();
    }
}

static inline int mandel(float c_re, float c_im, int count) {
    float z_re = c_re, z_im = c_im;
    int i;
    for (i = 0; i < count; ++i) {
        if (z_re * z_re + z_im * z_im > 4.f)
            break;

        float new_re = z_re * z_re - z_im * z_im;
        float new_im = 2.f * z_re * z_im;
        z_re = c_re + new_re;
        z_im = c_im + new_im;
    }
    return i;
}

