//
// Created by giovanni on 15/05/22.
//

#include "Kernel.h"

#define KERNEL_WIDTH 3

int kernelSize() {
    return KERNEL_WIDTH * KERNEL_WIDTH;
}

int kernelWidth() {
    return KERNEL_WIDTH;
}

float* createKernel(int type) {
    switch (type) {
        case kernelsType::boxBlur:
            return createBoxBlurKernel();
        case kernelsType::gaussianBlur:
            return createGaussianBlurKernel();
        case kernelsType::emboss:
            return createEmbossKernel();
        case kernelsType::outline:
            return createOutlineKernel();
        case kernelsType::sharpen:
            return createSharpenKernel();
        default:
            std::cerr << "Filter type not available." << std::endl;
            return nullptr;
    }
}
float* allocateEmptyKernel() {
    auto* kernel = (float*) malloc(sizeof(float) * kernelSize());
    return kernel;
}

float* createBoxBlurKernel() {
    auto* kernel = allocateEmptyKernel();
    for (int i = 0; i < kernelSize(); i++) {
        kernel[i] = 0.111111;
    }

    return kernel;
}

float* createEmbossKernel(){
    auto* kernel = allocateEmptyKernel();

    for (int i = 0; i < kernelSize(); i++) {
        if (i == 0)
            kernel[i] = -2.f;
        else if (i == kernelSize() - 1)
            kernel[i] = 2.f;
        else if (i == 2 || i == 6)
            kernel[i] = 0.f;
        else if (i == 1 || i == 3)
            kernel[i] = -1.f;
        else kernel[i] = 1;

    }

    return kernel;
}

float* createGaussianBlurKernel() {
    auto* kernel = allocateEmptyKernel();

    for (int i = 0; i < kernelSize(); ++i) {
        if (i == kernelSize()/2)
            kernel[i] = 1/4.f;
        else if (i % 2 == 0)
            kernel[i] = 1/16.f;
        else kernel[i] = 1/8.f;
    }

    return kernel;
}

float* createOutlineKernel() {
    auto* kernel = allocateEmptyKernel();

    for (int i = 0; i < kernelSize(); i++){
        if(i == kernelSize()/2)
            kernel[i] = 8.f;
        else
            kernel[i] = -1.f;
    }
    return kernel;
}

float* createSharpenKernel() {
    auto* kernel = allocateEmptyKernel();

    for (int i = 0; i < kernelSize(); ++i) {
        if (i == kernelSize()/2)
            kernel[i] = 5.f;
        else if (i % 2 == 0)
            kernel[i] = -1.f;
        else kernel[i] = 0.f;
    }

    return kernel;
}
