//
// Created by giovanni on 15/05/22.
//

#include "Convolution.h"
#include <iostream>

Image_t* convolution(Image_t* image, float* kernel, int kernelWidth){
    Image_t* processed = new_image(image->width, image->height, image->channels);

    float imagePixel;
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int kernelRadius = kernelWidth / 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;
    float* processedIter = processed->data;

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[(i * width + j) * channels + k] = 0;
                for (int y = -kernelRadius; y <= kernelRadius; y++) {
                    for (int x = -kernelRadius; x <= kernelRadius; x++) {
                        xOffset = j + x;
                        yOffset = i + y;
                        if (xOffset >= 0 && xOffset < width &&
                            yOffset >= 0 && yOffset < height) {
                            imagePixel = imageIter[(yOffset * width + xOffset) * channels + k];
                            kernelValue = kernelIter[(y + kernelRadius) * kernelWidth + x + kernelRadius];
                            processedIter[(i * width + j) * channels + k] += (imagePixel * kernelValue);
                        }
                    }
                }
            }
        }
    }
    return processed;
}

Image_t* convolutionNoBorder(Image_t* image, float* kernel, int kernelWidth){
    Image_t* processed = new_image(image->width - 2, image->height - 2, image->channels);

    float imagePixel;
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int kernelRadius = kernelWidth / 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;
    float* processedIter = processed->data;

    for (int i = 1; i < height-1; i++) {
        for (int j = 1; j < width-1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i - 1) * (width - 2) + (j - 1)) * channels + k] = 0;
                for (int y = -kernelRadius; y <= kernelRadius; y++) {
                    for (int x = -kernelRadius; x <= kernelRadius; x++) {
                        xOffset = j + x;
                        yOffset = i + y;
                        imagePixel = imageIter[(yOffset * width + xOffset) * channels + k];
                        kernelValue = kernelIter[(y + kernelRadius) * kernelWidth + x + kernelRadius];
                        processedIter[((i - 1) * (width - 2) + (j - 1)) * channels + k] += (imagePixel * kernelValue);

                    }
                }
            }
        }
    }
    return processed;
}

Image_t* convolutionUnrolling(Image_t* image, float* kernel, int kernelWidth){
    if(kernelWidth != 3) {
        std::cerr << "Kernel not supported, it only accepts 3x3 kernels" << std::endl;
        return nullptr;
    }

    Image_t* processed = new_image(image->width - 2, image->height - 2, image->channels);

    int width = image->width;
    int height = image->height;
    int channels = image->channels;

    float* imageIter = image->data;
    float* kernelIter = kernel;
    float* processedIter = processed->data;

    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * (width - 2) + (j-1)) * channels + k] =
                        imageIter[((i-1) * width + j-1) * channels + k] * kernelIter[0] +
                        imageIter[((i-1) * width + j) * channels + k] * kernelIter[1] +
                        imageIter[((i-1) * width + j+1) * channels + k] * kernelIter[2] +
                        imageIter[(i * width + j-1) * channels + k] * kernelIter[3] +
                        imageIter[(i * width + j) * channels + k] * kernelIter[4] +
                        imageIter[(i * width + j+1) * channels + k] * kernelIter[5] +
                        imageIter[((i+1) * width + j-1) * channels + k] * kernelIter[6] +
                        imageIter[((i+1) * width + j) * channels + k] * kernelIter[7] +
                        imageIter[((i+1) * width + j+1) * channels + k] * kernelIter[8];
            }
        }
    }

    return processed;
}

// it needs a processed image of full size, withBorders
Image_t* convolutionOMPNaive(Image_t* image, float* kernel, int kernelWidth, int nThreads){
    Image_t* processed = new_image(image->width, image->height, image->channels);

    float imagePixel;
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int kernelRadius = kernelWidth / 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;
    float* processedIter = processed->data;
#pragma omp parallel for default(none) \
    shared(width, height, channels, kernelRadius, kernelWidth, imageIter, kernelIter, processedIter) \
    firstprivate(xOffset, yOffset, imagePixel, kernelValue) \
    collapse(2) schedule(static) num_threads(nThreads)
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[(i * width + j) * channels + k] = 0;
                for (int y = -kernelRadius; y <= kernelRadius; y++) {
                    for (int x = -kernelRadius; x <= kernelRadius; x++) {
                        xOffset = j + x;
                        yOffset = i + y;
                        if (xOffset >= 0 && xOffset < width &&
                            yOffset >= 0 && yOffset < height) {
                            imagePixel = imageIter[(yOffset * width + xOffset) * channels + k];
                            kernelValue = kernelIter[(y + kernelRadius) * kernelWidth + x + kernelRadius];
                            processedIter[(i * width + j) * channels + k] += (imagePixel * kernelValue);
                        }
                    }
                }
            }
        }
    }

    return processed;
}

Image_t* convolutionOMPNaiveNoBorder(Image_t* image, float* kernel, int kernelWidth, int nThreads){
    Image_t* processed = new_image(image->width - 2, image->height - 2, image->channels);

    float imagePixel;
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int kernelRadius = kernelWidth / 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;
    float* processedIter = processed->data;
#pragma omp parallel for default(none) \
    shared(width, height, channels, kernelRadius, kernelWidth, imageIter, kernelIter, processedIter) \
    firstprivate(xOffset, yOffset, imagePixel, kernelValue) \
    collapse(2) schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * (width - 2) + (j-1)) * channels + k] = 0;
                for (int y = -kernelRadius; y <= kernelRadius; y++) {
                    for (int x = -kernelRadius; x <= kernelRadius; x++) {
                        xOffset = (j) + x;
                        yOffset = (i) + y;
                        imagePixel = imageIter[(yOffset * width + xOffset) * channels + k];
                        kernelValue = kernelIter[(y + kernelRadius) * kernelWidth + x + kernelRadius];
                        processedIter[((i-1) * (width - 2) + (j-1)) * channels + k] += (imagePixel * kernelValue);
                    }
                }
            }
        }
    }

    return processed;
}

// valgrind possibly lost leak, is it my problem, or it's openmp??
// I need to align memory access for simd??
// I need to use restrict on pointers for simd??
// If i use accum and then clamp, it is a read after write dependency, it can be used by simd??
// I could not do clamp here, but do it on the export PPM so i'll avoid a matrix scan only for doing the clamp (sum of convolution and assign it to the result matrix directly. or can i declare clamp as vectorialized? declare simd tipo
// Unrolling is done better if i use Structure of Arrays instead of array of structure like here
Image_t* convolutionOMPUnrollingSIMDWidth(Image_t* image, float* kernel, int kernelWidth, int nThreads){
    if(kernelWidth != 3) {
        std::cerr << "Kernel not supported, it only accepts 3x3 kernels" << std::endl;
        return nullptr;
    }

    Image_t* processed = new_image(image->width - 2, image->height - 2, image->channels);

    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int kernelRadius = kernelWidth / 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;
    float* processedIter = processed->data;

#pragma omp parallel for default(none) \
    shared(width, height, channels, kernelRadius, kernelWidth, imageIter, kernelIter, processedIter) \
    schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
#pragma omp simd
        for (int j = 1; j < width - 1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * (width - 2) + (j-1)) * channels + k] = imageIter[((i-1) * width + j-1) * channels + k] * kernelIter[0] +
                        imageIter[((i-1) * width + j) * channels + k] * kernelIter[1] +
                        imageIter[((i-1) * width + j+1) * channels + k] * kernelIter[2] +
                        imageIter[(i * width + j-1) * channels + k] * kernelIter[3] +
                        imageIter[(i * width + j) * channels + k] * kernelIter[4] +
                        imageIter[(i * width + j+1) * channels + k] * kernelIter[5] +
                        imageIter[((i+1) * width + j-1) * channels + k] * kernelIter[6] +
                        imageIter[((i+1) * width + j) * channels + k] * kernelIter[7] +
                        imageIter[((i+1) * width + j+1) * channels + k] * kernelIter[8];
            }
        }
    }

    return processed;
}

Image_t* convolutionOMPUnrollingSIMDChannels(Image_t* image, float* kernel, int kernelWidth, int nThreads){
    if(kernelWidth != 3) {
        std::cerr << "Kernel not supported, it only accepts 3x3 kernels" << std::endl;
        return nullptr;
    }

    Image_t* processed = new_image(image->width - 2, image->height - 2, image->channels);

    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int kernelRadius = kernelWidth / 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;
    float* processedIter = processed->data;

#pragma omp parallel for default(none) \
    shared(width, height, channels, kernelRadius, kernelWidth, imageIter, kernelIter, processedIter) \
    collapse(2) schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
#pragma omp simd
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * (width - 2) + (j-1)) * channels + k] = imageIter[((i-1) * width + j-1) * channels + k] * kernelIter[0] +
                        imageIter[((i-1) * width + j) * channels + k] * kernelIter[1] +
                        imageIter[((i-1) * width + j+1) * channels + k] * kernelIter[2] +
                        imageIter[(i * width + j-1) * channels + k] * kernelIter[3] +
                        imageIter[(i * width + j) * channels + k] * kernelIter[4] +
                        imageIter[(i * width + j+1) * channels + k] * kernelIter[5] +
                        imageIter[((i+1) * width + j-1) * channels + k] * kernelIter[6] +
                        imageIter[((i+1) * width + j) * channels + k] * kernelIter[7] +
                        imageIter[((i+1) * width + j+1) * channels + k] * kernelIter[8];
            }
        }
    }

    return processed;
}

Image_t* convolutionOMPUnrollingDoubleSIMD(Image_t* image, float* kernel, int kernelWidth, int nThreads){
    if(kernelWidth != 3) {
        std::cerr << "Kernel not supported, it only accepts 3x3 kernels" << std::endl;
        return nullptr;
    }

    Image_t* processed = new_image(image->width - 2, image->height - 2, image->channels);

    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int kernelRadius = kernelWidth / 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;
    float* processedIter = processed->data;

#pragma omp parallel for default(none) \
    shared(width, height, channels, kernelRadius, kernelWidth, imageIter, kernelIter, processedIter) \
    schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
#pragma omp simd
        for (int j = 1; j < width - 1; j++) {
#pragma omp ordered simd
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * (width - 2) + (j-1)) * channels + k] = imageIter[((i-1) * width + j-1) * channels + k] * kernelIter[0] +
                        imageIter[((i-1) * width + j) * channels + k] * kernelIter[1] +
                        imageIter[((i-1) * width + j+1) * channels + k] * kernelIter[2] +
                        imageIter[(i * width + j-1) * channels + k] * kernelIter[3] +
                        imageIter[(i * width + j) * channels + k] * kernelIter[4] +
                        imageIter[(i * width + j+1) * channels + k] * kernelIter[5] +
                        imageIter[((i+1) * width + j-1) * channels + k] * kernelIter[6] +
                        imageIter[((i+1) * width + j) * channels + k] * kernelIter[7] +
                        imageIter[((i+1) * width + j+1) * channels + k] * kernelIter[8];
            }
        }
    }

    return processed;
}

Image_t* convolutionOMPUnrollingSIMDCollapse(Image_t* image, float* kernel, int kernelWidth, int nThreads){
    if(kernelWidth != 3) {
        std::cerr << "Kernel not supported, it only accepts 3x3 kernels" << std::endl;
        return nullptr;
    }

    Image_t* processed = new_image(image->width - 2, image->height - 2, image->channels);

    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int kernelRadius = kernelWidth / 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;
    float* processedIter = processed->data;

#pragma omp parallel for default(none) \
    shared(width, height, channels, kernelRadius, kernelWidth, imageIter, kernelIter, processedIter) \
    schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
#pragma omp simd collapse(2)
        for (int j = 1; j < width - 1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * (width - 2) + (j-1)) * channels + k] = imageIter[((i-1) * width + j-1) * channels + k] * kernelIter[0] +
                        imageIter[((i-1) * width + j) * channels + k] * kernelIter[1] +
                        imageIter[((i-1) * width + j+1) * channels + k] * kernelIter[2] +
                        imageIter[(i * width + j-1) * channels + k] * kernelIter[3] +
                        imageIter[(i * width + j) * channels + k] * kernelIter[4] +
                        imageIter[(i * width + j+1) * channels + k] * kernelIter[5] +
                        imageIter[((i+1) * width + j-1) * channels + k] * kernelIter[6] +
                        imageIter[((i+1) * width + j) * channels + k] * kernelIter[7] +
                        imageIter[((i+1) * width + j+1) * channels + k] * kernelIter[8];
            }
        }
    }

    return processed;
}

Image_t* convolutionOMPUnrollingParallelForSIMD(Image_t* image, float* kernel, int kernelWidth, int nThreads){
    if(kernelWidth != 3) {
        std::cerr << "Kernel not supported, it only accepts 3x3 kernels" << std::endl;
        return nullptr;
    }

    Image_t* processed = new_image(image->width - 2, image->height - 2, image->channels);

    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int kernelRadius = kernelWidth / 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;
    float* processedIter = processed->data;

#pragma omp parallel for simd default(none) \
    shared(width, height, channels, kernelRadius, kernelWidth, imageIter, kernelIter, processedIter) \
    collapse(2) schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * (width - 2) + (j-1)) * channels + k] = imageIter[((i-1) * width + j-1) * channels + k] * kernelIter[0] +
                        imageIter[((i-1) * width + j) * channels + k] * kernelIter[1] +
                        imageIter[((i-1) * width + j+1) * channels + k] * kernelIter[2] +
                        imageIter[(i * width + j-1) * channels + k] * kernelIter[3] +
                        imageIter[(i * width + j) * channels + k] * kernelIter[4] +
                        imageIter[(i * width + j+1) * channels + k] * kernelIter[5] +
                        imageIter[((i+1) * width + j-1) * channels + k] * kernelIter[6] +
                        imageIter[((i+1) * width + j) * channels + k] * kernelIter[7] +
                        imageIter[((i+1) * width + j+1) * channels + k] * kernelIter[8];
            }
        }
    }

    return processed;
}

Image_t* convolutionOMPUnrollingSIMDWidthRestrict(Image_t* image, float* kernel, int kernelWidth, int nThreads){
    if(kernelWidth != 3) {
        std::cerr << "Kernel not supported, it only accepts 3x3 kernels" << std::endl;
        return nullptr;
    }

    Image_t* processed = new_image(image->width - 2, image->height - 2, image->channels);

    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int kernelRadius = kernelWidth / 2;

    float* __restrict imageIter = image->data;
    float* __restrict kernelIter = kernel;
    float* __restrict processedIter = processed->data;

#pragma omp parallel for default(none) \
    shared(width, height, channels, kernelRadius, kernelWidth, imageIter, kernelIter, processedIter) \
    schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
#pragma omp simd
        for (int j = 1; j < width - 1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * (width - 2) + (j-1)) * channels + k] =
                        imageIter[((i-1) * width + j-1) * channels + k] * kernelIter[0] +
                        imageIter[((i-1) * width + j) * channels + k] * kernelIter[1] +
                        imageIter[((i-1) * width + j+1) * channels + k] * kernelIter[2] +
                        imageIter[(i * width + j-1) * channels + k] * kernelIter[3] +
                        imageIter[(i * width + j) * channels + k] * kernelIter[4] +
                        imageIter[(i * width + j+1) * channels + k] * kernelIter[5] +
                        imageIter[((i+1) * width + j-1) * channels + k] * kernelIter[6] +
                        imageIter[((i+1) * width + j) * channels + k] * kernelIter[7] +
                        imageIter[((i+1) * width + j+1) * channels + k] * kernelIter[8];
            }
        }
    }

    return processed;
}