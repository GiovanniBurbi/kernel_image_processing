//
// Created by giovanni on 15/05/22.
//

#ifndef KERNEL_IMAGE_PROCESSING_CONVOLUTION_H
#define KERNEL_IMAGE_PROCESSING_CONVOLUTION_H

#include "image/Image.h"
#include "image/ImageSoA.h"

Image_t* convolution(Image_t* image, float* kernel);
ImageSoA_t* convolutionSoA(ImageSoA_t* image, float* kernel);
Image_t* convolutionUnrolling(Image_t* image, float* kernel);
Image_t* convolutionUnrollingKernel(Image_t* image, float* kernel);
Image_t* convolutionUnrollingChannels(Image_t* image, float* kernel);
ImageSoA_t* convolutionUnrollingSoA(ImageSoA_t* image, float* kernel);
Image_t* convolutionOMPNaive(Image_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPNaiveUnrollingChannels(Image_t* image, float* kernel, int nThreads);
ImageSoA_t* convolutionOMPNaiveSoA(ImageSoA_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrolling(Image_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrollingKernel(Image_t* image, float* kernel, int nThreads);
ImageSoA_t* convolutionOMPUnrollingSoA(ImageSoA_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrollingSIMDWidth(Image_t* image, float* kernel, int nThreads);
ImageSoA_t* convolutionOMPUnrollingSIMDWidthSoA(ImageSoA_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrollingSIMDChannels(Image_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrollingDoubleSIMD(Image_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrollingSIMDCollapse(Image_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrollingParallelForSIMD(Image_t* image, float* kernel, int nThreads);
Image_t* convolutionOMPUnrollingSIMDWidthRestrict(Image_t* image, float* kernel, int nThreads);
ImageSoA_t* convolutionOMPUnrollingSIMDWidthRestrictSoA(ImageSoA_t* image, float* kernel, int nThreads);

#endif //KERNEL_IMAGE_PROCESSING_CONVOLUTION_H
