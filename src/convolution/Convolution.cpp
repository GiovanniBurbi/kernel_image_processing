//
// Created by giovanni on 15/05/22.
//

#include "Convolution.h"
#include "utils/Utils.h"

void convolution(Image_t* image, Image_t* processed, float* kernel, int kernelWidth){
    float imagePixel;
    float kernelValue;
    float accum;
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
                accum = 0;
                for (int y = -kernelRadius; y <= kernelRadius; y++) {
                    for (int x = -kernelRadius; x <= kernelRadius; x++) {
                        xOffset = j + x;
                        yOffset = i + y;
                        if (xOffset >= 0 && xOffset < width &&
                            yOffset >= 0 && yOffset < height) {
                            imagePixel = imageIter[(yOffset * width + xOffset) * channels + k];
                            kernelValue = kernelIter[(y + kernelRadius) * kernelWidth + x + kernelRadius];
                            accum += (imagePixel * kernelValue);
                        }
                    }
                }
                processedIter[(i * width + j) * channels + k] = clamp(accum, 0, 1);
            }
        }
    }
}