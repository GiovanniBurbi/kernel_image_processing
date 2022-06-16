//
// Created by giovanni on 15/05/22.
//

#ifndef KERNEL_IMAGE_PROCESSING_KERNEL_H
#define KERNEL_IMAGE_PROCESSING_KERNEL_H

#include <iostream>

#define KERNEL_WIDTH 3
#define KERNEL_RADIUS KERNEL_WIDTH / 2

enum kernelsType {
    boxBlur = 0,
    gaussianBlur = 1,
    emboss = 2,
    outline = 3,
    sharpen = 4
};

float* createKernel(int type);
float* allocateEmptyKernel();
float* createBoxBlurKernel();
float* createEmbossKernel();
float* createGaussianBlurKernel();
float* createOutlineKernel();
float* createSharpenKernel();

#endif //KERNEL_IMAGE_PROCESSING_KERNEL_H
