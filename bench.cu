#include "cufft.h"
#include <iostream>

#define C2R 1
#define R2C 2
#define C2C 3
#define Z2D 5
#define D2Z 6
#define Z2Z 7
#define _FROMTO FROMTO

#if _FROMTO == Z2Z
#define TO_TYPE cufftDoubleComplex
#define FROM_TYPE cufftDoubleComplex
#define FROMTO_STR "double precision complex-to-complex"
#elif _FROMTO == D2Z
#define TO_TYPE cufftDoubleComplex
#define FROM_TYPE cufftDoubleReal
#define FROMTO_STR "double precision real-to-complex"
#elif _FROMTO == Z2D
#define TO_TYPE cufftDoubleReal
#define FROM_TYPE cufftDoubleComplex
#define FROMTO_STR "double precision complex-to-real"
#elif _FROMTO == C2C
#define TO_TYPE cufftComplex
#define FROM_TYPE cufftComplex
#define FROMTO_STR "single precision complex-to-complex"
#elif _FROMTO == R2C
#define TO_TYPE cufftComplex
#define FROM_TYPE cufftReal
#define FROMTO_STR "single precision real-to-complex"
#elif _FROMTO == C2R
#define TO_TYPE cufftReal
#define FROM_TYPE cufftComplex
#define FROMTO_STR "single precision complex-to-real"
#else
#error "FROMTO must be one of Z2Z, Z2D, D2Z, C2C, R2C and C2R"
#endif
template <class A, class B>
cufftResult_t CUFFTPLAN2D(cufftHandle *plan, int size_x, int size_y, A* in, B* out);

cufftResult_t CUFFTPLAN2D( cufftHandle *plan, int size_x, int size_y, 
                     cufftDoubleComplex* in, cufftDoubleComplex* out) {
      return cufftPlan2d(plan, size_x, size_y, CUFFT_Z2Z);
}

cufftResult_t CUFFTPLAN2D( cufftHandle *plan, int size_x, int size_y, 
                     cufftDoubleReal* in, cufftDoubleComplex* out) {
      return cufftPlan2d(plan, size_x, size_y, CUFFT_D2Z);
}

cufftResult_t CUFFTPLAN2D( cufftHandle *plan, int size_x, int size_y, 
                     cufftDoubleComplex* in, cufftDoubleReal* out) {
      return cufftPlan2d(plan, size_x, size_y, CUFFT_Z2D);
}

cufftResult_t CUFFTPLAN2D( cufftHandle *plan, int size_x, int size_y, 
                     cufftComplex* in, cufftComplex* out) {
      return cufftPlan2d(plan, size_x, size_y, CUFFT_C2C);
}

cufftResult_t CUFFTPLAN2D( cufftHandle *plan, int size_x, int size_y, 
                     cufftReal* in, cufftComplex* out) {
      return cufftPlan2d(plan, size_x, size_y, CUFFT_R2C);
}

cufftResult_t CUFFTPLAN2D( cufftHandle *plan, int size_x, int size_y, 
                     cufftComplex* in, cufftReal* out) {
      return cufftPlan2d(plan, size_x, size_y, CUFFT_C2R);
}

template <class A, class B>
cufftResult_t CUFFTEXEC(cufftHandle plan, A* in, B* out); 

cufftResult_t CUFFTEXEC (
                             cufftHandle plan, cufftDoubleComplex* in, cufftDoubleComplex* out) {
      return cufftExecZ2Z(plan, in, out, CUFFT_FORWARD);
}

cufftResult_t CUFFTEXEC(
                             cufftHandle plan, cufftDoubleReal* in, cufftDoubleComplex* out) {
      return cufftExecD2Z(plan, in, out);
}

cufftResult_t CUFFTEXEC(
                             cufftHandle plan, cufftDoubleComplex* in, cufftDoubleReal* out) {
      return cufftExecZ2D(plan, in, out);
}

cufftResult_t CUFFTEXEC(
                             cufftHandle plan, cufftComplex* in, cufftComplex* out) {
      return cufftExecC2C(plan, in, out, CUFFT_FORWARD);
}

cufftResult_t CUFFTEXEC(
                             cufftHandle plan, cufftReal* in, cufftComplex* out) {
      return cufftExecR2C(plan, in, out);
}

cufftResult_t CUFFTEXEC(
                             cufftHandle plan, cufftComplex* in, cufftReal* out) {
      return cufftExecC2R(plan, in, out);
}

int main(void) {

  int NX=10112, NY=10112;
  int size = NX*NY;
  float elapsed;
  cufftHandle plan;
  FROM_TYPE *data1;
  cudaMalloc(&data1, sizeof(FROM_TYPE)*NX*NY);
#ifndef INPLACE
  TO_TYPE *data2;
  cudaMalloc(&data2, sizeof(TO_TYPE)*NX*NY);
#endif

  cudaEvent_t start, stop;
  cudaEventCreate(&start); cudaEventCreate(&stop);

  cudaError_t err;
  err = cudaGetLastError();
  if (err) std::cout << "Error in initial copy" << std::endl;
  std::cin >> NX >> NY;
  std::cout << "**** " << FROMTO_STR << " ****" << std::endl;
  std::cout << "dx, dy, elapsed, Gcell/s, Gflps" << std::endl; 
#ifdef INPLACE 
#define TARGET data1
#else
#define TARGET data2
#endif
  while( NX != 0) {
    cufftResult_t r = CUFFTPLAN2D(&plan, NX, NY, data1, TARGET);
    cudaEventRecord(start);
    for (int z=0; z< 5; z++)
       if (!r) r = CUFFTEXEC(plan, data1, TARGET);
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsed, start, stop);
    cudaError_t err;
    err = cudaGetLastError();
    if (err) std::cout << NX << ", " << NY << " - Error " << err <<" : " <<
                        cudaGetErrorString(err) << std::endl;
    else if (r) std::cout << NX << ", " << NY << " - CUFFT Error " << r << 
                        std::endl;
    else std::cout << NX << ", " << NY << ", " << elapsed/5 << ", " 
              << 5*NX*NY/elapsed/1000/1000 << ", " << 5*5/elapsed/1000/1000*NX*NY*(log2(NX+0.000)+log2(NY+0.000)) << std::endl;
    cufftDestroy(plan);
    std::cin >> NX >> NY;
    if (NX*NY > size) {
       std::cout << "Reallocating to " << NX << " x " << NY << std::endl;
       cudaFree(data1); data1=0;
       cudaMalloc(&data1, sizeof(cufftDoubleComplex)*NX*NY);
       if(!data1) std::cout << "Failed to allocate data1!" << std::endl;
#ifndef INPLACE
       cudaFree(data2); data2=0;
       cudaMalloc(&data2, sizeof(cufftDoubleComplex)*NX*NY);
       if(!data2) std::cout << "Failed to allocate data2!\n" << std::endl;
#endif
       size = NX*NY;
    }
  }
  std::cout << "0, 0" << std::endl;
 
  //printf("(%d,%d) - Error %d: %s\n", NX, NY, err, cudaGetErrorString(err));

  cudaFree(data1); 
#ifndef INPLACE
  cudaFree(data2);
#endif
  cudaEventDestroy(start); cudaEventDestroy(stop);


  return 0;
}
