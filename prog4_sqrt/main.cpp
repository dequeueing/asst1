#include <stdio.h>
#include <algorithm>
#include <pthread.h>
#include <math.h>
#include <immintrin.h>  // AVX2 intrinsics

#include "CycleTimer.h"
#include "sqrt_ispc.h"

using namespace ispc;

extern void sqrtSerial(int N, float startGuess, float* values, float* output);

static void verifyResult(int N, float* result, float* gold) {
    for (int i=0; i<N; i++) {
        if (fabs(result[i] - gold[i]) > 1e-4) {
            printf("Error: [%d] Got %f expected %f\n", i, result[i], gold[i]);
        }
    }
}

// AVX2 implementation of the sqrt function
// Uses manual vectorization with explicit mask handling for control flow divergence
void sqrtAVX2(int N, float initialGuess, float* values, float* output) {
    const float kThreshold = 0.00001f;
    
    // AVX2 constants
    const __m256 threshold = _mm256_set1_ps(kThreshold);
    const __m256 initial_guess = _mm256_set1_ps(initialGuess);
    const __m256 three = _mm256_set1_ps(3.0f);
    const __m256 half = _mm256_set1_ps(0.5f);
    const __m256 one = _mm256_set1_ps(1.0f);
    
    int i = 0;
    
    // Process 8 floats at a time with AVX2
    for (i = 0; i <= N - 8; i += 8) {
        // Load 8 input values
        __m256 x = _mm256_loadu_ps(&values[i]);
        
        // Initialize guess for all 8 values
        __m256 guess = initial_guess;
        
        // Calculate initial predicate: |guess^2 * x - 1|
        __m256 guess_squared = _mm256_mul_ps(guess, guess);
        __m256 pred = _mm256_mul_ps(guess_squared, x);
        pred = _mm256_sub_ps(pred, one);
        // Absolute value using bit manipulation (clear sign bit)
        const __m256 abs_mask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF));
        pred = _mm256_and_ps(pred, abs_mask);
        
        // Create mask for values that need more iterations
        __m256 mask = _mm256_cmp_ps(pred, threshold, _CMP_GT_OQ);
        
        // Iterative Newton's method - continue while any lane needs more iterations
        while (!_mm256_testz_ps(mask, mask)) {  // while mask is not all zeros
            // Newton update: guess = (3*guess - x*guess^3) * 0.5
            __m256 guess_cubed = _mm256_mul_ps(guess_squared, guess);
            __m256 x_times_guess_cubed = _mm256_mul_ps(x, guess_cubed);
            __m256 three_times_guess = _mm256_mul_ps(three, guess);
            __m256 numerator = _mm256_sub_ps(three_times_guess, x_times_guess_cubed);
            __m256 new_guess = _mm256_mul_ps(numerator, half);
            
            // Update guess only for lanes that haven't converged
            guess = _mm256_blendv_ps(guess, new_guess, mask);
            
            // Recalculate predicate
            guess_squared = _mm256_mul_ps(guess, guess);
            pred = _mm256_mul_ps(guess_squared, x);
            pred = _mm256_sub_ps(pred, one);
            pred = _mm256_and_ps(pred, abs_mask);  // abs value
            
            // Update mask
            mask = _mm256_cmp_ps(pred, threshold, _CMP_GT_OQ);
        }
        
        // Final result: output = x * guess
        __m256 result = _mm256_mul_ps(x, guess);
        
        // Store result
        _mm256_storeu_ps(&output[i], result);
    }
    
    // Handle remaining elements with scalar code
    for (; i < N; i++) {
        float x = values[i];
        float guess = initialGuess;
        float pred = fabs(guess * guess * x - 1.0f);
        
        while (pred > kThreshold) {
            guess = (3.0f * guess - x * guess * guess * guess) * 0.5f;
            pred = fabs(guess * guess * x - 1.0f);
        }
        
        output[i] = x * guess;
    }
}

int main() {

    const unsigned int N = 20 * 1000 * 1000;
    const float initialGuess = 1.0f;

    float* values = new float[N];
    float* output = new float[N];
    float* gold = new float[N];

    for (unsigned int i=0; i<N; i++)
    {
        // TODO: CS149 students.  Attempt to change the values in the
        // array here to meet the instructions in the handout: we want
        // to you generate best and worse-case speedups
        
        // starter code populates array with random input values
        values[i] = .001f + 2.998f * static_cast<float>(rand()) / RAND_MAX;
    }

    // generate a gold version to check results
    for (unsigned int i=0; i<N; i++)
        gold[i] = sqrt(values[i]);

    //
    // And run the serial implementation 3 times, again reporting the
    // minimum time.
    //
    double minSerial = 1e30;
    for (int i = 0; i < 3; ++i) {
        double startTime = CycleTimer::currentSeconds();
        sqrtSerial(N, initialGuess, values, output);
        double endTime = CycleTimer::currentSeconds();
        minSerial = std::min(minSerial, endTime - startTime);
    }

    printf("[sqrt serial]:\t\t[%.3f] ms\n", minSerial * 1000);

    verifyResult(N, output, gold);

    //
    // Compute the image using the ispc implementation; report the minimum
    // time of three runs.
    //
    double minISPC = 1e30;
    for (int i = 0; i < 3; ++i) {
        double startTime = CycleTimer::currentSeconds();
        sqrt_ispc(N, initialGuess, values, output);
        double endTime = CycleTimer::currentSeconds();
        minISPC = std::min(minISPC, endTime - startTime);
    }

    printf("[sqrt ispc]:\t\t[%.3f] ms\n", minISPC * 1000);

    verifyResult(N, output, gold);

    // Clear out the buffer
    for (unsigned int i = 0; i < N; ++i)
        output[i] = 0;

    //
    // Tasking version of the ISPC code
    //
    double minTaskISPC = 1e30;
    for (int i = 0; i < 3; ++i) {
        double startTime = CycleTimer::currentSeconds();
        sqrt_ispc_withtasks(N, initialGuess, values, output);
        double endTime = CycleTimer::currentSeconds();
        minTaskISPC = std::min(minTaskISPC, endTime - startTime);
    }

    printf("[sqrt task ispc]:\t[%.3f] ms\n", minTaskISPC * 1000);

    verifyResult(N, output, gold);

    // Clear out the buffer
    for (unsigned int i = 0; i < N; ++i)
        output[i] = 0;

    //
    // Test our custom AVX2 implementation
    //
    double minAVX2 = 1e30;
    for (int i = 0; i < 3; ++i) {
        double startTime = CycleTimer::currentSeconds();
        sqrtAVX2(N, initialGuess, values, output);
        double endTime = CycleTimer::currentSeconds();
        minAVX2 = std::min(minAVX2, endTime - startTime);
    }

    printf("[sqrt AVX2]:\t\t[%.3f] ms\n", minAVX2 * 1000);

    verifyResult(N, output, gold);

    printf("\t\t\t\t(%.2fx speedup from ISPC)\n", minSerial/minISPC);
    printf("\t\t\t\t(%.2fx speedup from task ISPC)\n", minSerial/minTaskISPC);
    printf("\t\t\t\t(%.2fx speedup from AVX2)\n", minSerial/minAVX2);
    printf("\t\t\t\t(AVX2 vs ISPC: %.2fx)\n", minISPC/minAVX2);

    delete [] values;
    delete [] output;
    delete [] gold;

    return 0;
}
