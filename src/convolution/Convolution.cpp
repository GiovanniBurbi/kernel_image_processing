//
// Created by giovanni on 15/05/22.
//

#include "Convolution.h"

#define KERNEL_RADIUS 1
#define KERNEL_WIDTH 3


//Image_t* convolution(Image_t* image, float* kernel){
//    float imagePixelR;
//    float imagePixelG;
//    float imagePixelB;
//    float kernelValue;
//    int xOffset;
//    int yOffset;
//
//    int width = image->width;
//    int height = image->height;
//    int channels = image->channels;
//    int processedWidth = image->width - 2;
//
//    float* imageIter = image->data;
//    float* kernelIter = kernel;
//
//    Image_t* processed = new_image(processedWidth, image->height - 2, image->channels);
//    float* processedIter = processed->data;
//
//    for (int i = 1; i < height-1; i++) {
//        for (int j = 1; j < width-1; j++) {
//            processedIter[((i - 1) * processedWidth + (j - 1)) * channels ] = 0;
//            processedIter[((i - 1) * processedWidth + (j - 1)) * channels + 1] = 0;
//            processedIter[((i - 1) * processedWidth + (j - 1)) * channels + 2] = 0;
//            for (int y = -KERNEL_RADIUS; y <= KERNEL_RADIUS; y++) {
//                for (int x = -KERNEL_RADIUS; x <= KERNEL_RADIUS; x++) {
//                    xOffset = j + x;
//                    yOffset = i + y;
//                    imagePixelR = imageIter[(yOffset * width + xOffset) * channels];
//                    imagePixelG = imageIter[(yOffset * width + xOffset) * channels + 1];
//                    imagePixelB = imageIter[(yOffset * width + xOffset) * channels + 2];
//                    kernelValue = kernelIter[(y + KERNEL_RADIUS) * KERNEL_WIDTH + x + KERNEL_RADIUS];
//                    processedIter[((i - 1) * processedWidth + (j - 1)) * channels] += (imagePixelR * kernelValue);
//                    processedIter[((i - 1) * processedWidth + (j - 1)) * channels + 1] += (imagePixelG * kernelValue);
//                    processedIter[((i - 1) * processedWidth + (j - 1)) * channels + 2] += (imagePixelB * kernelValue);
//
//                }
//            }
//        }
//    }
//    return processed;
//}

 Image_t* convolution(Image_t* image, float* kernel){
    float imagePixel;
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int processedWidth = image->width - 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, image->height - 2, image->channels);
    float* processedIter = processed->data;

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

ImageSoA_t* convolutionSoA(ImageSoA_t* image, float* kernel){
    float imageRPixel;
    float imageGPixel;
    float imageBPixel;
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image->width;
    int height = image->height;
    int processedWidth = image->width - 2;

    float* imageRIter = image->r;
    float* imageGIter = image->g;
    float* imageBIter = image->b;
    float* kernelIter = kernel;

    ImageSoA_t* processed = new_imageSoA(processedWidth, image->height - 2);
    float* processedRIter = processed->r;
    float* processedGIter = processed->g;
    float* processedBIter = processed->b;

    for (int i = 1; i < height-1; i++) {
        for (int j = 1; j < width-1; j++) {
            processedRIter[(i-1) * processedWidth + (j-1)] = 0;
            processedGIter[(i-1) * processedWidth + (j-1)] = 0;
            processedBIter[(i-1) * processedWidth + (j-1)] = 0;
            for (int y = -KERNEL_RADIUS; y <= KERNEL_RADIUS; y++) {
                for (int x = -KERNEL_RADIUS; x <= KERNEL_RADIUS; x++) {
                    xOffset = j + x;
                    yOffset = i + y;
                    imageRPixel = imageRIter[(yOffset * width + xOffset)];
                    imageGPixel = imageGIter[(yOffset * width + xOffset)];
                    imageBPixel = imageBIter[(yOffset * width + xOffset)];
                    kernelValue = kernelIter[(y + KERNEL_RADIUS) * KERNEL_WIDTH + x + KERNEL_RADIUS];
                    processedRIter[((i - 1) * processedWidth + (j - 1))] += (imageRPixel * kernelValue);
                    processedGIter[((i - 1) * processedWidth + (j - 1))] += (imageGPixel * kernelValue);
                    processedBIter[((i - 1) * processedWidth + (j - 1))] += (imageBPixel * kernelValue);
                }
            }
        }
    }
    return processed;
}

Image_t* convolutionUnrolling(Image_t* image, float* kernel){
    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int processedWidth = image->width - 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, image->height - 2, image->channels);
    float* processedIter = processed->data;

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

// Image_t* convolutionUnrolling(Image_t* image, float* kernel){
//    int width = image->width;
//    int height = image->height;
//    int channels = image->channels;
//    int processedWidth = image->width - 2;
//
//    float* imageIter = image->data;
//    float* kernelIter = kernel;
//
//    Image_t* processed = new_image(processedWidth, image->height - 2, image->channels);
//    float* processedIter = processed->data;
//
//    for (int i = 1; i < height - 1; i++) {
//        for (int j = 1; j < width - 1; j++) {
//            for (int k = 0; k < channels; k++) {
//                processedIter[((i-1) * processedWidth + (j-1)) * channels + k] =
//                        imageIter[((i-1) * width + j-1) * channels + k] * kernelIter[0] +
//                        imageIter[((i-1) * width + j) * channels + k] * kernelIter[1] +
//                        imageIter[((i-1) * width + j+1) * channels + k] * kernelIter[2] +
//                        imageIter[(i * width + j-1) * channels + k] * kernelIter[3] +
//                        imageIter[(i * width + j) * channels + k] * kernelIter[4] +
//                        imageIter[(i * width + j+1) * channels + k] * kernelIter[5] +
//                        imageIter[((i+1) * width + j-1) * channels + k] * kernelIter[6] +
//                        imageIter[((i+1) * width + j) * channels + k] * kernelIter[7] +
//                        imageIter[((i+1) * width + j+1) * channels + k] * kernelIter[8];
//            }
//        }
//    }
//
//    return processed;
//}

ImageSoA_t* convolutionUnrollingSoA(ImageSoA_t* image, float* kernel){
    float imageRPixel;
    float imageGPixel;
    float imageBPixel;
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image->width;
    int height = image->height;
    int processedWidth = image->width - 2;

    float* imageRIter = image->r;
    float* imageGIter = image->g;
    float* imageBIter = image->b;
    float* kernelIter = kernel;

    ImageSoA_t* processed = new_imageSoA(processedWidth, image->height - 2);
    float* processedRIter = processed->r;
    float* processedGIter = processed->g;
    float* processedBIter = processed->b;

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

Image_t* convolutionOMPNaive(Image_t* image, float* kernel, int nThreads){
    float imagePixel;
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int processedWidth = image->width - 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, image->height - 2, image->channels);
    float* processedIter = processed->data;

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

Image_t* convolutionOMPNaiveUnrollingChannels(Image_t* image, float* kernel, int nThreads){
    float imagePixel1;
    float imagePixel2;
    float imagePixel3;
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int processedWidth = image->width - 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, image->height - 2, image->channels);
    float* processedIter = processed->data;

#pragma omp parallel for default(none) \
    shared(width, height, channels, imageIter, kernelIter, processedIter, processedWidth) \
    firstprivate(xOffset, yOffset, imagePixel1, imagePixel2, imagePixel3, kernelValue) \
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
                    imagePixel1 = imageIter[(yOffset * width + xOffset) * channels];
                    processedIter[((i-1) * processedWidth + (j-1)) * channels] += (imagePixel1 * kernelValue);
                    imagePixel2 = imageIter[(yOffset * width + xOffset) * channels + 1];
                    processedIter[((i-1) * processedWidth + (j-1)) * channels + 1] += (imagePixel2 * kernelValue);
                    imagePixel3 = imageIter[(yOffset * width + xOffset) * channels + 2];
                    processedIter[((i-1) * processedWidth + (j-1)) * channels + 2] += (imagePixel3 * kernelValue);
                }
            }
        }
    }

    return processed;
}

ImageSoA_t* convolutionOMPNaiveSoA(ImageSoA_t* image, float* kernel, int nThreads){
    float imageRPixel;
    float imageGPixel;
    float imageBPixel;
    float kernelValue;
    int xOffset;
    int yOffset;

    int width = image->width;
    int height = image->height;
    int processedWidth = image->width - 2;

    float* imageRIter = image->r;
    float* imageGIter = image->g;
    float* imageBIter = image->b;
    float* kernelIter = kernel;

    ImageSoA_t* processed = new_imageSoA(processedWidth, image->height - 2);
    float* processedRIter = processed->r;
    float* processedGIter = processed->g;
    float* processedBIter = processed->b;

#pragma omp parallel for default(none) \
    shared(width, height, imageRIter, imageGIter, imageBIter, \
    kernelIter, processedRIter, processedGIter, processedBIter, processedWidth) \
    firstprivate(xOffset, yOffset, imageRPixel, imageGPixel, imageBPixel, \
    kernelValue) \
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
                    imageRPixel = imageRIter[(yOffset * width + xOffset)];
                    imageGPixel = imageGIter[(yOffset * width + xOffset)];
                    imageBPixel = imageBIter[(yOffset * width + xOffset)];
                    kernelValue = kernelIter[(y + KERNEL_RADIUS) * KERNEL_WIDTH + x + KERNEL_RADIUS];
                    processedRIter[((i-1) * processedWidth + (j-1))] += (imageRPixel * kernelValue);
                    processedBIter[((i-1) * processedWidth + (j-1))] += (imageGPixel * kernelValue);
                    processedGIter[((i-1) * processedWidth + (j-1))] += (imageBPixel * kernelValue);
                }
            }
        }
    }

    return processed;
}

Image_t* convolutionOMPUnrollingSIMDWidth(Image_t* image, float* kernel, int nThreads){
    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int processedWidth = image->width - 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, image->height - 2, image->channels);
    float* processedIter = processed->data;

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

ImageSoA_t* convolutionOMPUnrollingSIMDWidthSoA(ImageSoA_t* image, float* kernel, int nThreads){
    int width = image->width;
    int height = image->height;
    int processedWidth = image->width - 2;

    float* imageRIter = image->r;
    float* imageGIter = image->g;
    float* imageBIter = image->b;
    float* kernelIter = kernel;

    ImageSoA_t* processed = new_imageSoA(processedWidth, image->height - 2);
    float* processedRIter = processed->r;
    float* processedGIter = processed->g;
    float* processedBIter = processed->b;

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

Image_t* convolutionOMPUnrolling(Image_t* image, float* kernel, int nThreads){
    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int processedWidth = image->width - 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, image->height - 2, image->channels);
    float* processedIter = processed->data;

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

ImageSoA_t* convolutionOMPUnrollingSoA(ImageSoA_t* image, float* kernel, int nThreads){
    int width = image->width;
    int height = image->height;
    int processedWidth = image->width - 2;

    float* imageRIter = image->r;
    float* imageGIter = image->g;
    float* imageBIter = image->b;
    float* kernelIter = kernel;

    ImageSoA_t* processed = new_imageSoA(processedWidth, image->height - 2);
    float* processedRIter = processed->r;
    float* processedGIter = processed->g;
    float* processedBIter = processed->b;

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

Image_t* convolutionOMPUnrollingSIMDChannels(Image_t* image, float* kernel, int nThreads){
    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int processedWidth = image->width - 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, image->height - 2, image->channels);
    float* processedIter = processed->data;

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

Image_t* convolutionOMPUnrollingDoubleSIMD(Image_t* image, float* kernel, int nThreads){
    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int processedWidth = image->width - 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, image->height - 2, image->channels);
    float* processedIter = processed->data;

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

Image_t* convolutionOMPUnrollingSIMDCollapse(Image_t* image, float* kernel, int nThreads){
    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int processedWidth = image->width - 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, image->height - 2, image->channels);
    float* processedIter = processed->data;

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

Image_t* convolutionOMPUnrollingParallelForSIMD(Image_t* image, float* kernel, int nThreads){
    int width = image->width;
    int height = image->height;
    int channels = image->channels;
    int processedWidth = image->width - 2;

    float* imageIter = image->data;
    float* kernelIter = kernel;

    Image_t* processed = new_image(processedWidth, image->height - 2, image->channels);
    float* processedIter = processed->data;

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

    ImageSoA_t* processed = new_imageSoA(processedWidth, image->height - 2);
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