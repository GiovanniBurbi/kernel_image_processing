//
// Created by giovanni on 15/05/22.
//

#include "Convolution.h"
#include <iostream>

#define KERNEL_RADIUS 1
#define KERNEL_WIDTH 3
#define PIXEL_LOST_PER_AXIS 2
#define RGB_CHANNELS 3

/*
 * Sequential convolution naive implementation
 * It works for images with any number of channels and any kernels
 * */
 Image_t* convolution(Image_t* image, float* kernel){
    float imagePixel;
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int channels = image_getChannels(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    float* imageIter = image_getData(image);
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, height - PIXEL_LOST_PER_AXIS, channels);
    float* processedIter = image_getData(processed);

    for (int i = 1; i < height-1; i++) {
        for (int j = 1; j < width-1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i - 1) * processedWidth + (j - 1)) * channels + k] = 0;
                for (int y = -KERNEL_RADIUS; y <= KERNEL_RADIUS; y++) {
                    for (int x = -KERNEL_RADIUS; x <= KERNEL_RADIUS; x++) {
                        xOffset = j + x;
                        yOffset = i + y;
                        imagePixel = imageIter[(yOffset * width + xOffset) * channels + k];
                        kernelValue = kernelIter[(y + KERNEL_RADIUS) * KERNEL_WIDTH + x + KERNEL_RADIUS];
                        processedIter[((i - 1) * processedWidth + (j - 1)) * channels + k] += (imagePixel * kernelValue);

                    }
                }
            }
        }
    }
    return processed;
}

/*
 * Sequential convolution implementation with SoA data layout
 * It works for images with three channels and any kernel
 * */
ImageSoA_t* convolutionSoA(ImageSoA_t* image, float* kernel){
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    if (image_getChannels(image) != RGB_CHANNELS) {
        std::cerr << "This convolution implementation works only with images with three channels " << std::endl;
        return nullptr;
    }

    float* imageRIter = image_getR(image);
    float* imageGIter = image_getG(image);
    float* imageBIter = image_getB(image);
    float* kernelIter = kernel;

    ImageSoA_t* processed = new_imageSoA(processedWidth, height - PIXEL_LOST_PER_AXIS, RGB_CHANNELS);
    float* processedRIter = image_getR(processed);
    float* processedGIter = image_getG(processed);
    float* processedBIter = image_getB(processed);

    for (int i = 1; i < height-1; i++) {
        for (int j = 1; j < width-1; j++) {
            processedRIter[(i-1) * processedWidth + (j-1)] = 0;
            processedGIter[(i-1) * processedWidth + (j-1)] = 0;
            processedBIter[(i-1) * processedWidth + (j-1)] = 0;

            for (int y = -KERNEL_RADIUS; y <= KERNEL_RADIUS; y++) {
                for (int x = -KERNEL_RADIUS; x <= KERNEL_RADIUS; x++) {
                    xOffset = j + x;
                    yOffset = i + y;

                    kernelValue = kernelIter[(y + KERNEL_RADIUS) * KERNEL_WIDTH + x + KERNEL_RADIUS];

                    processedRIter[((i - 1) * processedWidth + (j - 1))] += (imageRIter[(yOffset * width + xOffset)] * kernelValue);
                    processedGIter[((i - 1) * processedWidth + (j - 1))] += (imageGIter[(yOffset * width + xOffset)] * kernelValue);
                    processedBIter[((i - 1) * processedWidth + (j - 1))] += (imageBIter[(yOffset * width + xOffset)] * kernelValue);
                }
            }
        }
    }
    return processed;
}

/*
 * Sequential convolution implementation with loop unrolling
 * It works for images with three channels and 3x3 kernels
 * */
Image_t* convolutionUnrolling(Image_t* image, float* kernel){
    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int channels = image_getChannels(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    if (channels != 3) {
        std::cerr << "This convolution implementation works only with images with three channels " << std::endl;
        return nullptr;
    }

    float* imageIter = image_getData(image);
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, height - PIXEL_LOST_PER_AXIS, RGB_CHANNELS);
    float* processedIter = image_getData(processed);

    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            processedIter[((i-1) * processedWidth + (j-1)) * channels] =
                    imageIter[((i-1) * width + j-1) * channels] * kernelIter[0] +
                    imageIter[((i-1) * width + j) * channels] * kernelIter[1] +
                    imageIter[((i-1) * width + j+1) * channels] * kernelIter[2] +
                    imageIter[(i * width + j-1) * channels] * kernelIter[3] +
                    imageIter[(i * width + j) * channels] * kernelIter[4] +
                    imageIter[(i * width + j+1) * channels] * kernelIter[5] +
                    imageIter[((i+1) * width + j-1) * channels] * kernelIter[6] +
                    imageIter[((i+1) * width + j) * channels] * kernelIter[7] +
                    imageIter[((i+1) * width + j+1) * channels] * kernelIter[8];

            processedIter[((i-1) * processedWidth + (j-1)) * channels + 1] =
                    imageIter[((i-1) * width + j-1) * channels + 1] * kernelIter[0] +
                    imageIter[((i-1) * width + j) * channels + 1] * kernelIter[1] +
                    imageIter[((i-1) * width + j+1) * channels + 1] * kernelIter[2] +
                    imageIter[(i * width + j-1) * channels + 1] * kernelIter[3] +
                    imageIter[(i * width + j) * channels + 1] * kernelIter[4] +
                    imageIter[(i * width + j+1) * channels + 1] * kernelIter[5] +
                    imageIter[((i+1) * width + j-1) * channels + 1] * kernelIter[6] +
                    imageIter[((i+1) * width + j) * channels + 1] * kernelIter[7] +
                    imageIter[((i+1) * width + j+1) * channels + 1] * kernelIter[8];

            processedIter[((i-1) * processedWidth + (j-1)) * channels + 2] =
                    imageIter[((i-1) * width + j-1) * channels + 2] * kernelIter[0] +
                    imageIter[((i-1) * width + j) * channels + 2] * kernelIter[1] +
                    imageIter[((i-1) * width + j+1) * channels + 2] * kernelIter[2] +
                    imageIter[(i * width + j-1) * channels + 2] * kernelIter[3] +
                    imageIter[(i * width + j) * channels + 2] * kernelIter[4] +
                    imageIter[(i * width + j+1) * channels + 2] * kernelIter[5] +
                    imageIter[((i+1) * width + j-1) * channels + 2] * kernelIter[6] +
                    imageIter[((i+1) * width + j) * channels + 2] * kernelIter[7] +
                    imageIter[((i+1) * width + j+1) * channels + 2] * kernelIter[8];
        }
    }

    return processed;
}

/*
 * Sequential convolution implementation with loop unrolling of kernel only
 * It works for images with any number of channels and 3x3 kernels
 * */
 Image_t* convolutionUnrollingKernel(Image_t* image, float* kernel){
    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int channels = image_getChannels(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    float* imageIter = image_getData(image);
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, height - PIXEL_LOST_PER_AXIS, channels);
    float* processedIter = image_getData(processed);

    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * processedWidth + (j-1)) * channels + k] =
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

/*
 * Sequential convolution implementation with loop unrolling of channels only
 * It works only for images with three channels and any kernel
 * */
Image_t* convolutionUnrollingChannels(Image_t* image, float* kernel){
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int channels = image_getChannels(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    if (channels != RGB_CHANNELS) {
        std::cerr << "This convolution implementation works only with images with three channels " << std::endl;
        return nullptr;
    }

    float* imageIter = image_getData(image);
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, height - PIXEL_LOST_PER_AXIS, RGB_CHANNELS);
    float* processedIter = image_getData(processed);

    for (int i = 1; i < height-1; i++) {
        for (int j = 1; j < width-1; j++) {
            processedIter[((i - 1) * processedWidth + (j - 1)) * channels ] = 0;
            processedIter[((i - 1) * processedWidth + (j - 1)) * channels + 1] = 0;
            processedIter[((i - 1) * processedWidth + (j - 1)) * channels + 2] = 0;

            for (int y = -KERNEL_RADIUS; y <= KERNEL_RADIUS; y++) {
                for (int x = -KERNEL_RADIUS; x <= KERNEL_RADIUS; x++) {
                    xOffset = j + x;
                    yOffset = i + y;

                    kernelValue = kernelIter[(y + KERNEL_RADIUS) * KERNEL_WIDTH + x + KERNEL_RADIUS];
                    processedIter[((i - 1) * processedWidth + (j - 1)) * channels] += (imageIter[(yOffset * width + xOffset) * channels] * kernelValue);
                    processedIter[((i - 1) * processedWidth + (j - 1)) * channels + 1] += (imageIter[(yOffset * width + xOffset) * channels + 1] * kernelValue);
                    processedIter[((i - 1) * processedWidth + (j - 1)) * channels + 2] += (imageIter[(yOffset * width + xOffset) * channels + 2] * kernelValue);

                }
            }
        }
    }
    return processed;
}

/*
 * Sequential convolution implementation with loop unrolling and SoA data layout
 * It works only for images with three channels and 3x3 kernels
 * */
ImageSoA_t* convolutionUnrollingSoA(ImageSoA_t* image, float* kernel){
    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    if (image_getChannels(image) != RGB_CHANNELS) {
        std::cerr << "This convolution implementation works only with images with three channels " << std::endl;
        return nullptr;
    }

    float* imageRIter = image_getR(image);
    float* imageGIter = image_getG(image);
    float* imageBIter = image_getB(image);
    float* kernelIter = kernel;

    ImageSoA_t* processed = new_imageSoA(processedWidth, height - PIXEL_LOST_PER_AXIS, RGB_CHANNELS);
    float* processedRIter = image_getR(processed);
    float* processedGIter = image_getG(processed);
    float* processedBIter = image_getB(processed);

    for (int i = 1; i < height-1; i++) {
        for (int j = 1; j < width-1; j++) {
            processedRIter[((i-1) * processedWidth + (j-1))] =
                    imageRIter[((i-1) * width + j-1)] * kernelIter[0] +
                    imageRIter[((i-1) * width + j)] * kernelIter[1] +
                    imageRIter[((i-1) * width + j+1)] * kernelIter[2] +
                    imageRIter[(i * width + j-1)] * kernelIter[3] +
                    imageRIter[(i * width + j)] * kernelIter[4] +
                    imageRIter[(i * width + j+1)] * kernelIter[5] +
                    imageRIter[((i+1) * width + j-1)] * kernelIter[6] +
                    imageRIter[((i+1) * width + j)] * kernelIter[7] +
                    imageRIter[((i+1) * width + j+1)] * kernelIter[8];

            processedGIter[((i-1) * processedWidth + (j-1))] =
                    imageGIter[((i-1) * width + j-1)] * kernelIter[0] +
                    imageGIter[((i-1) * width + j)] * kernelIter[1] +
                    imageGIter[((i-1) * width + j+1)] * kernelIter[2] +
                    imageGIter[(i * width + j-1)] * kernelIter[3] +
                    imageGIter[(i * width + j)] * kernelIter[4] +
                    imageGIter[(i * width + j+1)] * kernelIter[5] +
                    imageGIter[((i+1) * width + j-1)] * kernelIter[6] +
                    imageGIter[((i+1) * width + j)] * kernelIter[7] +
                    imageGIter[((i+1) * width + j+1)] * kernelIter[8];

            processedBIter[((i-1) * processedWidth + (j-1))] =
                    imageBIter[((i-1) * width + j-1)] * kernelIter[0] +
                    imageBIter[((i-1) * width + j)] * kernelIter[1] +
                    imageBIter[((i-1) * width + j+1)] * kernelIter[2] +
                    imageBIter[(i * width + j-1)] * kernelIter[3] +
                    imageBIter[(i * width + j)] * kernelIter[4] +
                    imageBIter[(i * width + j+1)] * kernelIter[5] +
                    imageBIter[((i+1) * width + j-1)] * kernelIter[6] +
                    imageBIter[((i+1) * width + j)] * kernelIter[7] +
                    imageBIter[((i+1) * width + j+1)] * kernelIter[8];
        }
    }
    return processed;
}

/*
 * Parallel convolution naive implementation with openMP
 * It works for images with any number of channels and any kernels
 * */
Image_t* convolutionOMPNaive(Image_t* image, float* kernel, int nThreads){
    float imagePixel;
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int channels = image_getChannels(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    float* imageIter = image_getData(image);
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, height - PIXEL_LOST_PER_AXIS, channels);
    float* processedIter = image_getData(processed);

#pragma omp parallel for default(none) \
    shared(width, height, channels, imageIter, kernelIter, processedIter, processedWidth) \
    firstprivate(xOffset, yOffset, imagePixel, kernelValue) \
    collapse(2) schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * processedWidth + (j-1)) * channels + k] = 0;
                for (int y = -KERNEL_RADIUS; y <= KERNEL_RADIUS; y++) {
                    for (int x = -KERNEL_RADIUS; x <= KERNEL_RADIUS; x++) {
                        xOffset = (j) + x;
                        yOffset = (i) + y;
                        imagePixel = imageIter[(yOffset * width + xOffset) * channels + k];
                        kernelValue = kernelIter[(y + KERNEL_RADIUS) * KERNEL_WIDTH + x + KERNEL_RADIUS];
                        processedIter[((i-1) * processedWidth + (j-1)) * channels + k] += (imagePixel * kernelValue);
                    }
                }
            }
        }
    }

    return processed;
}


/*
 * Parallel convolution implementation with openMP and loop unrolling channels only
 * It works only for images with three channels and any kernels
 * */
Image_t* convolutionOMPUnrollingChannels(Image_t* image, float* kernel, int nThreads){
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int channels = image_getChannels(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    if (channels != RGB_CHANNELS) {
        std::cerr << "This convolution implementation works only with images with three channels " << std::endl;
        return nullptr;
    }

    float* imageIter = image_getData(image);
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, height - PIXEL_LOST_PER_AXIS, RGB_CHANNELS);
    float* processedIter = image_getData(processed);

#pragma omp parallel for default(none) \
    shared(width, height, channels, imageIter, kernelIter, processedIter, processedWidth) \
    firstprivate(xOffset, yOffset, kernelValue) \
    collapse(2) schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            processedIter[((i-1) * processedWidth + (j-1)) * channels] = 0;
            processedIter[((i-1) * processedWidth + (j-1)) * channels + 1] = 0;
            processedIter[((i-1) * processedWidth + (j-1)) * channels + 2] = 0;

            for (int y = -KERNEL_RADIUS; y <= KERNEL_RADIUS; y++) {
                for (int x = -KERNEL_RADIUS; x <= KERNEL_RADIUS; x++) {
                    xOffset = j + x;
                    yOffset = i + y;
                    kernelValue = kernelIter[(y + KERNEL_RADIUS) * KERNEL_WIDTH + x + KERNEL_RADIUS];

                    processedIter[((i-1) * processedWidth + (j-1)) * channels] += (imageIter[(yOffset * width + xOffset) * channels] * kernelValue);
                    processedIter[((i-1) * processedWidth + (j-1)) * channels + 1] += (imageIter[(yOffset * width + xOffset) * channels + 1] * kernelValue);
                    processedIter[((i-1) * processedWidth + (j-1)) * channels + 2] += (imageIter[(yOffset * width + xOffset) * channels + 2] * kernelValue);
                }
            }
        }
    }

    return processed;
}


/*
 * Parallel convolution implementation with openMP and SoA data layout
 * It works only for images with three channels and any kernels
 * */
ImageSoA_t* convolutionOMPNaiveSoA(ImageSoA_t* image, float* kernel, int nThreads){
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    if (image_getChannels(image) != RGB_CHANNELS) {
        std::cerr << "This convolution implementation works only with images with three channels " << std::endl;
        return nullptr;
    }

    float* imageRIter = image_getR(image);
    float* imageGIter = image_getG(image);
    float* imageBIter = image_getB(image);
    float* kernelIter = kernel;

    ImageSoA_t* processed = new_imageSoA(processedWidth, height - PIXEL_LOST_PER_AXIS, RGB_CHANNELS);
    float* processedRIter = image_getR(processed);
    float* processedGIter = image_getG(processed);
    float* processedBIter = image_getB(processed);

#pragma omp parallel for default(none) \
    shared(width, height, imageRIter, imageGIter, imageBIter, \
    kernelIter, processedRIter, processedGIter, processedBIter, processedWidth) \
    firstprivate(xOffset, yOffset, kernelValue) \
    collapse(2) schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            processedRIter[((i-1) * processedWidth + (j-1))] = 0;
            processedGIter[((i-1) * processedWidth + (j-1))] = 0;
            processedBIter[((i-1) * processedWidth + (j-1))] = 0;

            for (int y = -KERNEL_RADIUS; y <= KERNEL_RADIUS; y++) {
                for (int x = -KERNEL_RADIUS; x <= KERNEL_RADIUS; x++) {
                    xOffset = (j) + x;
                    yOffset = (i) + y;
                    kernelValue = kernelIter[(y + KERNEL_RADIUS) * KERNEL_WIDTH + x + KERNEL_RADIUS];

                    processedRIter[((i-1) * processedWidth + (j-1))] += (imageRIter[(yOffset * width + xOffset)] * kernelValue);
                    processedBIter[((i-1) * processedWidth + (j-1))] += (imageGIter[(yOffset * width + xOffset)] * kernelValue);
                    processedGIter[((i-1) * processedWidth + (j-1))] += (imageBIter[(yOffset * width + xOffset)] * kernelValue);
                }
            }
        }
    }

    return processed;
}

/*
 * Parallel convolution implementation with openMP and loop unrolling kernel only
 * It works for images with any number of channels and 3x3 kernels
 * */
Image_t* convolutionOMPUnrollingSIMDWidth(Image_t* image, float* kernel, int nThreads){
    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int channels = image_getChannels(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    float* imageIter = image_getData(image);
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, height - PIXEL_LOST_PER_AXIS, channels);
    float* processedIter = image_getData(processed);

#pragma omp parallel for default(none) \
    shared(width, height, channels, imageIter, kernelIter, processedIter, processedWidth) \
    schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
#pragma omp simd
        for (int j = 1; j < width - 1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * processedWidth + (j-1)) * channels + k] = imageIter[((i-1) * width + j-1) * channels + k] * kernelIter[0] +
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

/*
 * Parallel convolution implementation with openMP, loop unrolling and SoA data layout
 * It works for images with three channels and 3x3 kernels
 * */
ImageSoA_t* convolutionOMPUnrollingSIMDWidthSoA(ImageSoA_t* image, float* kernel, int nThreads){
    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    if (image_getChannels(image) != RGB_CHANNELS) {
        std::cerr << "This convolution implementation works only with images with three channels " << std::endl;
        return nullptr;
    }

    float* imageRIter = image_getR(image);
    float* imageGIter = image_getG(image);
    float* imageBIter = image_getB(image);
    float* kernelIter = kernel;

    ImageSoA_t* processed = new_imageSoA(processedWidth, height - PIXEL_LOST_PER_AXIS, RGB_CHANNELS);
    float* processedRIter = image_getR(processed);
    float* processedGIter = image_getG(processed);
    float* processedBIter = image_getB(processed);

#pragma omp parallel for default(none) \
    shared(width, height, imageRIter, imageGIter, imageBIter, \
    kernelIter, processedRIter, processedGIter, processedBIter, processedWidth) \
    schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
#pragma omp simd
        for (int j = 1; j < width - 1; j++) {
            processedRIter[((i-1) * processedWidth + (j-1))] =
                    imageRIter[((i-1) * width + j-1)] * kernelIter[0] +
                    imageRIter[((i-1) * width + j)] * kernelIter[1] +
                    imageRIter[((i-1) * width + j+1)] * kernelIter[2] +
                    imageRIter[(i * width + j-1)] * kernelIter[3] +
                    imageRIter[(i * width + j)] * kernelIter[4] +
                    imageRIter[(i * width + j+1)] * kernelIter[5] +
                    imageRIter[((i+1) * width + j-1)] * kernelIter[6] +
                    imageRIter[((i+1) * width + j)] * kernelIter[7] +
                    imageRIter[((i+1) * width + j+1)] * kernelIter[8];

            processedGIter[((i-1) * processedWidth + (j-1))] =
                    imageGIter[((i-1) * width + j-1)] * kernelIter[0] +
                    imageGIter[((i-1) * width + j)] * kernelIter[1] +
                    imageGIter[((i-1) * width + j+1)] * kernelIter[2] +
                    imageGIter[(i * width + j-1)] * kernelIter[3] +
                    imageGIter[(i * width + j)] * kernelIter[4] +
                    imageGIter[(i * width + j+1)] * kernelIter[5] +
                    imageGIter[((i+1) * width + j-1)] * kernelIter[6] +
                    imageGIter[((i+1) * width + j)] * kernelIter[7] +
                    imageGIter[((i+1) * width + j+1)] * kernelIter[8];

            processedBIter[((i-1) * processedWidth + (j-1))] =
                    imageBIter[((i-1) * width + j-1)] * kernelIter[0] +
                    imageBIter[((i-1) * width + j)] * kernelIter[1] +
                    imageBIter[((i-1) * width + j+1)] * kernelIter[2] +
                    imageBIter[(i * width + j-1)] * kernelIter[3] +
                    imageBIter[(i * width + j)] * kernelIter[4] +
                    imageBIter[(i * width + j+1)] * kernelIter[5] +
                    imageBIter[((i+1) * width + j-1)] * kernelIter[6] +
                    imageBIter[((i+1) * width + j)] * kernelIter[7] +
                    imageBIter[((i+1) * width + j+1)] * kernelIter[8];
        }
    }

    return processed;
}

/*
 * Parallel convolution implementation with openMP and loop unrolling
 * It works for images with three channels and 3x3 kernels
 * */
Image_t* convolutionOMPUnrolling(Image_t* image, float* kernel, int nThreads){
    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int channels = image_getChannels(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    if (channels != RGB_CHANNELS) {
        std::cerr << "This convolution implementation works only with images with three channels " << std::endl;
        return nullptr;
    }

    float* imageIter = image_getData(image);
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, height - PIXEL_LOST_PER_AXIS, RGB_CHANNELS);
    float* processedIter = image_getData(processed);

#pragma omp parallel for default(none) \
    shared(width, height, channels, imageIter, kernelIter, processedIter, processedWidth) \
    collapse(2) schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            processedIter[((i-1) * processedWidth + (j-1)) * channels] =
                    imageIter[((i-1) * width + j-1) * channels] * kernelIter[0] +
                    imageIter[((i-1) * width + j) * channels] * kernelIter[1] +
                    imageIter[((i-1) * width + j+1) * channels] * kernelIter[2] +
                    imageIter[(i * width + j-1) * channels] * kernelIter[3] +
                    imageIter[(i * width + j) * channels] * kernelIter[4] +
                    imageIter[(i * width + j+1) * channels] * kernelIter[5] +
                    imageIter[((i+1) * width + j-1) * channels] * kernelIter[6] +
                    imageIter[((i+1) * width + j) * channels] * kernelIter[7] +
                    imageIter[((i+1) * width + j+1) * channels] * kernelIter[8];

            processedIter[((i-1) * processedWidth + (j-1)) * channels + 1] =
                    imageIter[((i-1) * width + j-1) * channels + 1] * kernelIter[0] +
                    imageIter[((i-1) * width + j) * channels + 1] * kernelIter[1] +
                    imageIter[((i-1) * width + j+1) * channels + 1] * kernelIter[2] +
                    imageIter[(i * width + j-1) * channels + 1] * kernelIter[3] +
                    imageIter[(i * width + j) * channels + 1] * kernelIter[4] +
                    imageIter[(i * width + j+1) * channels + 1] * kernelIter[5] +
                    imageIter[((i+1) * width + j-1) * channels + 1] * kernelIter[6] +
                    imageIter[((i+1) * width + j) * channels + 1] * kernelIter[7] +
                    imageIter[((i+1) * width + j+1) * channels + 1] * kernelIter[8];

            processedIter[((i-1) * processedWidth + (j-1)) * channels + 2] =
                    imageIter[((i-1) * width + j-1) * channels + 2] * kernelIter[0] +
                    imageIter[((i-1) * width + j) * channels + 2] * kernelIter[1] +
                    imageIter[((i-1) * width + j+1) * channels + 2] * kernelIter[2] +
                    imageIter[(i * width + j-1) * channels + 2] * kernelIter[3] +
                    imageIter[(i * width + j) * channels + 2] * kernelIter[4] +
                    imageIter[(i * width + j+1) * channels + 2] * kernelIter[5] +
                    imageIter[((i+1) * width + j-1) * channels + 2] * kernelIter[6] +
                    imageIter[((i+1) * width + j) * channels + 2] * kernelIter[7] +
                    imageIter[((i+1) * width + j+1) * channels + 2] * kernelIter[8];
        }
    }

    return processed;
}


/*
 * Parallel convolution implementation with openMP and loop unrolling kernel only
 * It works for images with any number of channels and 3x3 kernels
 * */
Image_t* convolutionOMPUnrollingKernel(Image_t* image, float* kernel, int nThreads){
    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int channels = image_getChannels(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    float* imageIter = image_getData(image);
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, height - PIXEL_LOST_PER_AXIS, channels);
    float* processedIter = image_getData(processed);

#pragma omp parallel for default(none) \
    shared(width, height, channels, imageIter, kernelIter, processedIter, processedWidth) \
    collapse(2) schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * processedWidth + (j-1)) * channels + k] =
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

/*
 * Parallel convolution implementation with openMP, loop unrolling and SoA data layout
 * It works for images with three channels and 3x3 kernels
 * */
ImageSoA_t* convolutionOMPUnrollingSoA(ImageSoA_t* image, float* kernel, int nThreads){
    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    if (image_getChannels(image) != RGB_CHANNELS) {
        std::cerr << "This convolution implementation works only with images with three channels " << std::endl;
        return nullptr;
    }

    float* imageRIter = image_getR(image);
    float* imageGIter = image_getG(image);
    float* imageBIter = image_getB(image);
    float* kernelIter = kernel;

    ImageSoA_t* processed = new_imageSoA(processedWidth, height - PIXEL_LOST_PER_AXIS, RGB_CHANNELS);
    float* processedRIter = image_getR(processed);
    float* processedGIter = image_getG(processed);
    float* processedBIter = image_getB(processed);

#pragma omp parallel for default(none) \
    shared(width, height, imageRIter, imageGIter, imageBIter, \
    kernelIter, processedRIter, processedGIter, processedBIter, processedWidth) \
    collapse(2) schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            processedRIter[((i-1) * processedWidth + (j-1))] =
                    imageRIter[((i-1) * width + j-1)] * kernelIter[0] +
                    imageRIter[((i-1) * width + j)] * kernelIter[1] +
                    imageRIter[((i-1) * width + j+1)] * kernelIter[2] +
                    imageRIter[(i * width + j-1)] * kernelIter[3] +
                    imageRIter[(i * width + j)] * kernelIter[4] +
                    imageRIter[(i * width + j+1)] * kernelIter[5] +
                    imageRIter[((i+1) * width + j-1)] * kernelIter[6] +
                    imageRIter[((i+1) * width + j)] * kernelIter[7] +
                    imageRIter[((i+1) * width + j+1)] * kernelIter[8];

            processedGIter[((i-1) * processedWidth + (j-1))] =
                    imageGIter[((i-1) * width + j-1)] * kernelIter[0] +
                    imageGIter[((i-1) * width + j)] * kernelIter[1] +
                    imageGIter[((i-1) * width + j+1)] * kernelIter[2] +
                    imageGIter[(i * width + j-1)] * kernelIter[3] +
                    imageGIter[(i * width + j)] * kernelIter[4] +
                    imageGIter[(i * width + j+1)] * kernelIter[5] +
                    imageGIter[((i+1) * width + j-1)] * kernelIter[6] +
                    imageGIter[((i+1) * width + j)] * kernelIter[7] +
                    imageGIter[((i+1) * width + j+1)] * kernelIter[8];

            processedBIter[((i-1) * processedWidth + (j-1))] =
                    imageBIter[((i-1) * width + j-1)] * kernelIter[0] +
                    imageBIter[((i-1) * width + j)] * kernelIter[1] +
                    imageBIter[((i-1) * width + j+1)] * kernelIter[2] +
                    imageBIter[(i * width + j-1)] * kernelIter[3] +
                    imageBIter[(i * width + j)] * kernelIter[4] +
                    imageBIter[(i * width + j+1)] * kernelIter[5] +
                    imageBIter[((i+1) * width + j-1)] * kernelIter[6] +
                    imageBIter[((i+1) * width + j)] * kernelIter[7] +
                    imageBIter[((i+1) * width + j+1)] * kernelIter[8];
        }
    }

    return processed;
}

/*
 * Parallel convolution implementation with openMP and loop unrolling kernel only
 * It works for images with any number of channels and 3x3 kernels
 * */
Image_t* convolutionOMPUnrollingSIMDChannels(Image_t* image, float* kernel, int nThreads){
    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int channels = image_getChannels(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    float* imageIter = image_getData(image);
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, height - PIXEL_LOST_PER_AXIS, channels);
    float* processedIter = image_getData(processed);

#pragma omp parallel for default(none) \
    shared(width, height, channels, imageIter, kernelIter, processedIter, processedWidth) \
    collapse(2) schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
#pragma omp simd
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * processedWidth + (j-1)) * channels + k] = imageIter[((i-1) * width + j-1) * channels + k] * kernelIter[0] +
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

/*
 * Parallel convolution implementation with openMP and loop unrolling kernel only
 * It works for images with any number of channels and 3x3 kernels
 * */
Image_t* convolutionOMPUnrollingDoubleSIMD(Image_t* image, float* kernel, int nThreads){
    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int channels = image_getChannels(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    float* imageIter = image_getData(image);
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, height - PIXEL_LOST_PER_AXIS, channels);
    float* processedIter = image_getData(processed);

#pragma omp parallel for default(none) \
    shared(width, height, channels, imageIter, kernelIter, processedIter, processedWidth) \
    schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
#pragma omp simd
        for (int j = 1; j < width - 1; j++) {
#pragma omp ordered simd
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * processedWidth + (j-1)) * channels + k] = imageIter[((i-1) * width + j-1) * channels + k] * kernelIter[0] +
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

/*
 * Parallel convolution implementation with openMP and loop unrolling kernel only
 * It works for images with any number of channels and 3x3 kernels
 * */
Image_t* convolutionOMPUnrollingSIMDCollapse(Image_t* image, float* kernel, int nThreads){
    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int channels = image_getChannels(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    float* imageIter = image_getData(image);
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, height - PIXEL_LOST_PER_AXIS, channels);
    float* processedIter = image_getData(processed);

#pragma omp parallel for default(none) \
    shared(width, height, channels, imageIter, kernelIter, processedIter, processedWidth) \
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

/*
 * Parallel convolution implementation with openMP and loop unrolling kernel only
 * It works for images with any number of channels and 3x3 kernels
 * */
Image_t* convolutionOMPUnrollingParallelForSIMD(Image_t* image, float* kernel, int nThreads){
    int width = image_getWidth(image);
    int height = image_getHeight(image);
    int channels = image_getChannels(image);
    int processedWidth = width - PIXEL_LOST_PER_AXIS;

    float* imageIter = image_getData(image);
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, height - PIXEL_LOST_PER_AXIS, channels);
    float* processedIter = image_getData(processed);

#pragma omp parallel for simd default(none) \
    shared(width, height, channels, imageIter, kernelIter, processedIter, processedWidth) \
    collapse(2) schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
        for (int j = 1; j < width - 1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * processedWidth + (j-1)) * channels + k] = imageIter[((i-1) * width + j-1) * channels + k] * kernelIter[0] +
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

Image_t* convolutionOMPUnrollingSIMDWidthRestrict(Image_t* image, float* kernel, int nThreads){
    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int processedWidth = image->width - 2;

    float* __restrict imageIter = image->data;
    float* __restrict kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, image->height - 2, image->channels);
    float* __restrict processedIter = processed->data;

#pragma omp parallel for default(none) \
    shared(width, height, channels, imageIter, kernelIter, processedIter, processedWidth) \
    schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
#pragma omp simd
        for (int j = 1; j < width - 1; j++) {
            for (int k = 0; k < channels; k++) {
                processedIter[((i-1) * processedWidth + (j-1)) * channels + k] =
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

ImageSoA_t* convolutionOMPUnrollingSIMDWidthRestrictSoA(ImageSoA_t* image, float* kernel, int nThreads){
    int width = image->width;
    int height = image->height;
    int processedWidth = image->width - 2;

    float* __restrict imageRIter = image->r;
    float* __restrict imageGIter = image->g;
    float* __restrict imageBIter = image->b;
    float* __restrict kernelIter = kernel;

    ImageSoA_t* processed = new_imageSoA(processedWidth, image->height - 2, 3);
    float* __restrict processedRIter = processed->r;
    float* __restrict processedGIter = processed->g;
    float* __restrict processedBIter = processed->b;

#pragma omp parallel for default(none) \
    shared(width, height, imageRIter, imageGIter, imageBIter, \
    kernelIter, processedRIter, processedGIter, processedBIter, processedWidth) \
    schedule(static) num_threads(nThreads)
    for (int i = 1; i < height - 1; i++) {
#pragma omp simd
        for (int j = 1; j < width - 1; j++) {
            processedRIter[((i-1) * processedWidth + (j-1))] =
                    imageRIter[((i-1) * width + j-1)] * kernelIter[0] +
                    imageRIter[((i-1) * width + j)] * kernelIter[1] +
                    imageRIter[((i-1) * width + j+1)] * kernelIter[2] +
                    imageRIter[(i * width + j-1)] * kernelIter[3] +
                    imageRIter[(i * width + j)] * kernelIter[4] +
                    imageRIter[(i * width + j+1)] * kernelIter[5] +
                    imageRIter[((i+1) * width + j-1)] * kernelIter[6] +
                    imageRIter[((i+1) * width + j)] * kernelIter[7] +
                    imageRIter[((i+1) * width + j+1)] * kernelIter[8];

            processedGIter[((i-1) * processedWidth + (j-1))] =
                    imageGIter[((i-1) * width + j-1)] * kernelIter[0] +
                    imageGIter[((i-1) * width + j)] * kernelIter[1] +
                    imageGIter[((i-1) * width + j+1)] * kernelIter[2] +
                    imageGIter[(i * width + j-1)] * kernelIter[3] +
                    imageGIter[(i * width + j)] * kernelIter[4] +
                    imageGIter[(i * width + j+1)] * kernelIter[5] +
                    imageGIter[((i+1) * width + j-1)] * kernelIter[6] +
                    imageGIter[((i+1) * width + j)] * kernelIter[7] +
                    imageGIter[((i+1) * width + j+1)] * kernelIter[8];

            processedBIter[((i-1) * processedWidth + (j-1))] =
                    imageBIter[((i-1) * width + j-1)] * kernelIter[0] +
                    imageBIter[((i-1) * width + j)] * kernelIter[1] +
                    imageBIter[((i-1) * width + j+1)] * kernelIter[2] +
                    imageBIter[(i * width + j-1)] * kernelIter[3] +
                    imageBIter[(i * width + j)] * kernelIter[4] +
                    imageBIter[(i * width + j+1)] * kernelIter[5] +
                    imageBIter[((i+1) * width + j-1)] * kernelIter[6] +
                    imageBIter[((i+1) * width + j)] * kernelIter[7] +
                    imageBIter[((i+1) * width + j+1)] * kernelIter[8];
        }
    }

    return processed;
}