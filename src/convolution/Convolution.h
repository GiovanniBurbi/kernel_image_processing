//
// Created by giovanni on 15/05/22.
//

#ifndef KERNEL_IMAGE_PROCESSING_CONVOLUTION_H
#define KERNEL_IMAGE_PROCESSING_CONVOLUTION_H

#include "image/Image.h"
#include "image/ImageSoA.h"

// Sequential Naive
Image_t* convolution(Image_t* image, const float* __restrict__ kernel);
// Sequential Naive SoA
ImageSoA_t* convolutionSoA(ImageSoA_t* image, const float* __restrict__ kernel);
// Sequential Unrolling kernel and channels
Image_t* convolutionUnrolling(Image_t* image, const float* __restrict__ kernel);
// Sequential Unrolling kernel only
Image_t* convolutionUnrollingKernel(Image_t* image, const float* __restrict__ kernel);
// Sequential Unrolling channels only
Image_t* convolutionUnrollingChannels(Image_t* image, const float* __restrict__ kernel);
// Sequential Unrolling kernel and channels SoA
ImageSoA_t* convolutionUnrollingSoA(ImageSoA_t* image, const float* __restrict__ kernel);

// Parallel with openMP naive
Image_t* convolutionOMPNaive(Image_t* image, const float* __restrict__ kernel, int nThreads);
// Parallel with openMP naive SoA (like unrolling channels)
ImageSoA_t* convolutionOMPNaiveSoA(ImageSoA_t* image, const float* __restrict__ kernel, int nThreads);
// Parallel with openMP naive Unrolling only channels
Image_t* convolutionOMPUnrollingChannels(Image_t* image, const float* __restrict__ kernel, int nThreads);
// Parallel with openMP naive Unrolling only kernel
Image_t* convolutionOMPUnrollingKernel(Image_t* image, const float* __restrict__ kernel, int nThreads);
// Parallel with openMP naive Unrolling kernel and channels
Image_t* convolutionOMPUnrolling(Image_t* image, const float* __restrict__ kernel, int nThreads);
// Parallel with openMP naive Unrolling kernel and channels SoA
ImageSoA_t* convolutionOMPUnrollingSoA(ImageSoA_t* image, const float* __restrict__ kernel, int nThreads);

// Parallel with openMP Unrolling kernel and simd width
Image_t* convolutionOMPUnrollingSIMDWidth(Image_t* image, const float* __restrict__ kernel, int nThreads);
// Parallel with openMP Unrolling kernel and simd width SoA (unrolling channels)
ImageSoA_t* convolutionOMPUnrollingSIMDWidthSoA(ImageSoA_t* image, const float* __restrict__ kernel, int nThreads);

// Parallel with openMP Unrolling kernel and simd channels
Image_t* convolutionOMPUnrollingSIMDChannels(Image_t* image, const float* __restrict__ kernel, int nThreads);
// Parallel with openMP Unrolling kernel double simd width and channels
Image_t* convolutionOMPUnrollingDoubleSIMD(Image_t* image, const float* __restrict__ kernel, int nThreads);
// Parallel with openMP Unrolling kernel simd collapse width and channels
Image_t* convolutionOMPUnrollingSIMDCollapse(Image_t* image, const float* __restrict__ kernel, int nThreads);
// Parallel with openMP Unrolling kernel simd inside parallel directive
Image_t* convolutionOMPUnrollingParallelForSIMD(Image_t* image, const float* __restrict__ kernel, int nThreads);

#endif //KERNEL_IMAGE_PROCESSING_CONVOLUTION_H
