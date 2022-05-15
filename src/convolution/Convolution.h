//
// Created by giovanni on 15/05/22.
//

#ifndef KERNEL_IMAGE_PROCESSING_CONVOLUTION_H
#define KERNEL_IMAGE_PROCESSING_CONVOLUTION_H

#include "image/Image.h"

void convolution(Image_t* image, Image_t* processed, float* kernel, int kernelWidth);
#endif //KERNEL_IMAGE_PROCESSING_CONVOLUTION_H
