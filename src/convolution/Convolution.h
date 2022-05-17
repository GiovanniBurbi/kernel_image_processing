//
// Created by giovanni on 15/05/22.
//

#ifndef KERNEL_IMAGE_PROCESSING_CONVOLUTION_H
#define KERNEL_IMAGE_PROCESSING_CONVOLUTION_H

#include "image/Image.h"

Image_t* convolution(Image_t* image, float* kernel);
Image_t* convolutionUnrolling(Image_t* image, float* kernel);
Image_t* convolutionOMPNaive(Image_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrolling(Image_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrollingSIMDWidth(Image_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrollingSIMDChannels(Image_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrollingDoubleSIMD(Image_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrollingSIMDCollapse(Image_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrollingParallelForSIMD(Image_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrollingSIMDWidthRestrict(Image_t* image, float* kernel, int nThreads);

#endif //KERNEL_IMAGE_PROCESSING_CONVOLUTION_H
