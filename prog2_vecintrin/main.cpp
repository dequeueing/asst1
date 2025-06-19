#include <stdio.h>
#include <algorithm>
#include <getopt.h>
#include <math.h>
#include "CS149intrin.h"
#include "logger.h"
using namespace std;

#define EXP_MAX 10

Logger CS149Logger;

void usage(const char* progname);
void initValue(float* values, int* exponents, float* output, float* gold, unsigned int N);
void absSerial(float* values, float* output, int N);
void absVector(float* values, float* output, int N);
void clampedExpSerial(float* values, int* exponents, float* output, int N);
void clampedExpVector(float* values, int* exponents, float* output, int N);
float arraySumSerial(float* values, int N);
float arraySumVector(float* values, int N);
bool verifyResult(float* values, int* exponents, float* output, float* gold, int N);

int wow(int argc, char * argv[]) {
  // try absVector
  
  float values[] = {1.0f, -2.0f, 3.0f, -4.0f, 5.0f, -6.0f, 7.0f, -8.0f, -9.0f};
  float output[10];

  absVector(values, output, 9);
  printf("Abs Vector Result:\n");
  for (int i = 0; i < 10; i++) {
    printf("%f ", output[i]);
  }
  printf("\n");
  return 0;
}

int main(int argc, char * argv[]) {
  int N = 16;
  bool printLog = false;

  // parse commandline options ////////////////////////////////////////////
  int opt;
  static struct option long_options[] = {
    {"size", 1, 0, 's'},
    {"log", 0, 0, 'l'},
    {"help", 0, 0, '?'},
    {0 ,0, 0, 0}
  };

  while ((opt = getopt_long(argc, argv, "s:l?", long_options, NULL)) != EOF) {

    switch (opt) {
      case 's':
        N = atoi(optarg);
        if (N <= 0) {
          printf("Error: Workload size is set to %d (<0).\n", N);
          return -1;
        }
        break;
      case 'l':
        printLog = true;
        break;
      case '?':
      default:
        usage(argv[0]);
        return 1;
    }
  }


  float* values = new float[N+VECTOR_WIDTH];
  int* exponents = new int[N+VECTOR_WIDTH];
  float* output = new float[N+VECTOR_WIDTH];
  float* gold = new float[N+VECTOR_WIDTH];
  initValue(values, exponents, output, gold, N);

  clampedExpSerial(values, exponents, gold, N);
  clampedExpVector(values, exponents, output, N);

  //absSerial(values, gold, N);
  //absVector(values, output, N);

  printf("\e[1;31mCLAMPED EXPONENT\e[0m (required) \n");
  bool clampedCorrect = verifyResult(values, exponents, output, gold, N);
  if (printLog) CS149Logger.printLog();
  CS149Logger.printStats();

  printf("************************ Result Verification *************************\n");
  if (!clampedCorrect) {
    printf("@@@ Failed!!!\n");
  } else {
    printf("Passed!!!\n");
  }

  printf("\n\e[1;31mARRAY SUM\e[0m (bonus) \n");
  if (N % VECTOR_WIDTH == 0) {
    float sumGold = arraySumSerial(values, N);
    float sumOutput = arraySumVector(values, N);
    float epsilon = 0.1;
    bool sumCorrect = abs(sumGold - sumOutput) < epsilon * 2;
    if (!sumCorrect) {
      printf("Expected %f, got %f\n.", sumGold, sumOutput);
      printf("@@@ Failed!!!\n");
    } else {
      printf("Passed!!!\n");
    }
  } else {
    printf("Must have N %% VECTOR_WIDTH == 0 for this problem (VECTOR_WIDTH is %d)\n", VECTOR_WIDTH);
  }

  delete [] values;
  delete [] exponents;
  delete [] output;
  delete [] gold;

  return 0;
}

void usage(const char* progname) {
  printf("Usage: %s [options]\n", progname);
  printf("Program Options:\n");
  printf("  -s  --size <N>     Use workload size N (Default = 16)\n");
  printf("  -l  --log          Print vector unit execution log\n");
  printf("  -?  --help         This message\n");
}

void initValue(float* values, int* exponents, float* output, float* gold, unsigned int N) {

  for (unsigned int i=0; i<N+VECTOR_WIDTH; i++)
  {
    // random input values
    values[i] = -1.f + 4.f * static_cast<float>(rand()) / RAND_MAX;
    exponents[i] = rand() % EXP_MAX;
    output[i] = 0.f;
    gold[i] = 0.f;
  }

}

bool verifyResult(float* values, int* exponents, float* output, float* gold, int N) {
  int incorrect = -1;
  float epsilon = 0.00001;
  for (int i=0; i<N+VECTOR_WIDTH; i++) {
    if ( abs(output[i] - gold[i]) > epsilon ) {
      incorrect = i;
      break;
    }
  }

  if (incorrect != -1) {
    if (incorrect >= N)
      printf("You have written to out of bound value!\n");
    printf("Wrong calculation at value[%d]!\n", incorrect);
    printf("value  = ");
    for (int i=0; i<N; i++) {
      printf("% f ", values[i]);
    } printf("\n");

    printf("exp    = ");
    for (int i=0; i<N; i++) {
      printf("% 9d ", exponents[i]);
    } printf("\n");

    printf("output = ");
    for (int i=0; i<N; i++) {
      printf("% f ", output[i]);
    } printf("\n");

    printf("gold   = ");
    for (int i=0; i<N; i++) {
      printf("% f ", gold[i]);
    } printf("\n");
    return false;
  }
  printf("Results matched with answer!\n");
  return true;
}

// computes the absolute value of all elements in the input array
// values, stores result in output
void absSerial(float* values, float* output, int N) {
  for (int i=0; i<N; i++) {
    float x = values[i];
    if (x < 0) {
      output[i] = -x;
    } else {
      output[i] = x;
    }
  }
}


// implementation of absSerial() above, but it is vectorized using CS149 intrinsics
void absVector(float* values, float* output, int N) {
  __cs149_vec_float x;
  __cs149_vec_float result;
  __cs149_vec_float zero = _cs149_vset_float(0.f);
  __cs149_mask maskAll, maskIsNegative, maskIsNotNegative;

//  Note: Take a careful look at this loop indexing.  This example
//  code is not guaranteed to work when (N % VECTOR_WIDTH) != 0.
//  Why is that the case?
  for (int i=0; i<N; i+=VECTOR_WIDTH) {

    // All ones
    maskAll = _cs149_init_ones();

    // All zeros
    maskIsNegative = _cs149_init_ones(0);

    // Load vector of values from contiguous memory addresses
    _cs149_vload_float(x, values+i, maskAll);               // x = values[i];

    // Set mask according to predicate
    _cs149_vlt_float(maskIsNegative, x, zero, maskAll);     // if (x < 0) {

    // Execute instruction using mask ("if" clause)
    _cs149_vsub_float(result, zero, x, maskIsNegative);      //   output[i] = -x;

    // Inverse maskIsNegative to generate "else" mask
    maskIsNotNegative = _cs149_mask_not(maskIsNegative);     // } else {

    // Execute instruction ("else" clause)
    _cs149_vload_float(result, values+i, maskIsNotNegative); //   output[i] = x; }

    // Write results back to memory
    _cs149_vstore_float(output+i, result, maskAll);
  }
}


// accepts an array of values and an array of exponents
//
// For each element, compute values[i]^exponents[i] and clamp value to
// 9.999.  Store result in output.
void clampedExpSerial(float* values, int* exponents, float* output, int N) {
  for (int i=0; i<N; i++) {
    float x = values[i];
    int y = exponents[i];
    if (y == 0) {
      output[i] = 1.f;
    } else {
      float result = x;
      int count = y - 1;
      while (count > 0) {
        result *= x;
        count--;
      }
      if (result > 9.999999f) {
        result = 9.999999f;
      }
      output[i] = result;
    }
  }
}

void clampedExpVector(float* values, int* exponents, float* output, int N) {

  //
  // CS149 STUDENTS TODO: Implement your vectorized version of
  // clampedExpSerial() here.
  //
  // Your solution should work for any value of
  // N and VECTOR_WIDTH, not just when VECTOR_WIDTH divides N
  //

  // Declare vector variables
  __cs149_vec_float x, result, clamp_val, one_val;
  __cs149_vec_int y, count, zero_int, one_int;
  __cs149_mask maskAll, maskExpZero, maskNotExpZero, maskCountPos, maskToClamp;
  
  // Initialize constant vectors
  clamp_val = _cs149_vset_float(9.999999f);
  one_val = _cs149_vset_float(1.0f);
  zero_int = _cs149_vset_int(0);
  one_int = _cs149_vset_int(1);

  // iterate each width of vectors and convert values into vector
  for (int i = 0; i < N; i += VECTOR_WIDTH) {
    
    // Handle case where N is not divisible by VECTOR_WIDTH
    int remaining = N - i;
    if (remaining >= VECTOR_WIDTH) {
      maskAll = _cs149_init_ones();
    } else {
      maskAll = _cs149_init_ones(remaining);
    }
    
    // Load input vectors
    _cs149_vload_float(x, values + i, maskAll);
    _cs149_vload_int(y, exponents + i, maskAll);
    
    // handle special case of exponent == 0
    _cs149_veq_int(maskExpZero, y, zero_int, maskAll);
    maskNotExpZero = _cs149_mask_not(maskExpZero);
    maskNotExpZero = _cs149_mask_and(maskNotExpZero, maskAll);
    
    // Set result to 1.0 for exponent == 0 cases
    _cs149_vset_float(result, 1.0f, maskExpZero);
    
    // vec1: the value, vec2: corresponding exponent - 1
    // For exponent != 0 cases, initialize result = x and count = y - 1
    _cs149_vmove_float(result, x, maskNotExpZero);
    _cs149_vsub_int(count, y, one_int, maskNotExpZero);
    
    // continue if there is exponent != 0
    // for each iteration, do the multiplication if the mask is not zero
    while (_cs149_cntbits(maskNotExpZero) > 0) {
      // Check which lanes still have count > 0
      _cs149_vgt_int(maskCountPos, count, zero_int, maskNotExpZero);
      
      // If no lanes have count > 0, break
      if (_cs149_cntbits(maskCountPos) == 0) {
        break;
      }
      
      // Multiply result by x for lanes with count > 0
      _cs149_vmult_float(result, result, x, maskCountPos);
      
      // Decrement count for lanes with count > 0
      _cs149_vsub_int(count, count, one_int, maskCountPos);
    }
    
    // Clamp result to 9.999999 if exceeded
    _cs149_vgt_float(maskToClamp, result, clamp_val, maskAll);
    _cs149_vmove_float(result, clamp_val, maskToClamp);
    
    // Store results
    _cs149_vstore_float(output + i, result, maskAll);
  }
}

// returns the sum of all elements in values
float arraySumSerial(float* values, int N) {
  float sum = 0;
  for (int i=0; i<N; i++) {
    sum += values[i];
  }

  return sum;
}

// returns the sum of all elements in values
// You can assume N is a multiple of VECTOR_WIDTH
// You can assume VECTOR_WIDTH is a power of 2
float arraySumVector(float* values, int N) {
  
  //
  // CS149 STUDENTS TODO: Implement your vectorized version of arraySumSerial here
  //
  
  __cs149_vec_float sum_vec, temp_vec;
  __cs149_mask maskAll = _cs149_init_ones();
  
  // Initialize sum vector to zeros
  sum_vec = _cs149_vset_float(0.0f);
  
  // Step 1: Sum all elements into the vector (N/VECTOR_WIDTH operations)
  for (int i = 0; i < N; i += VECTOR_WIDTH) {
    _cs149_vload_float(temp_vec, values + i, maskAll);
    _cs149_vadd_float(sum_vec, sum_vec, temp_vec, maskAll);
  }
  
  // Step 2: Tree reduction using hadd and interleave operations
  // We need log2(VECTOR_WIDTH) steps to reduce to single value
  
  __cs149_vec_float reduced_vec = sum_vec;
  
  // First hadd: pairs adjacent elements
  _cs149_hadd_float(reduced_vec, reduced_vec);
  
  // Continue reduction for any VECTOR_WIDTH >= 4
  // Each iteration reduces the problem size by half
  for (int step_width = VECTOR_WIDTH / 2; step_width >= 2; step_width /= 2) {
    // Interleave to rearrange elements for next reduction step
    _cs149_interleave_float(reduced_vec, reduced_vec);
    // Add adjacent pairs again
    _cs149_hadd_float(reduced_vec, reduced_vec);
  }
  
  // The final sum should now be in reduced_vec.value[0]
  return reduced_vec.value[0];
}

// 1,2,3,4,5,6,7,8
// 3,3,7,7,11,11,15,15
// 3,7,11,15,3,7,11,15
// 10,10,26,26,10,10,26,26
// 10,26,10,26,10,26,10,26
// 36,36,36,36,36,36,36,36