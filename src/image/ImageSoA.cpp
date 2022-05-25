//
// Created by giovanni on 18/05/22.
//

#include "ImageSoA.h"
#include <iostream>

ImageSoA_t* new_ImageSoA1(int width, int height, float *r, float *g, float *b) {
    ImageSoA_t* img;

    img = (ImageSoA_t*) malloc(sizeof(ImageSoA_t));

    image_setWidth(img, width);
    image_setHeight(img, height);

    image_setR(img, r);
    image_setG(img, g);
    image_setB(img, b);
    return img;
}

ImageSoA_t* new_imageSoA(int width, int height) {
    float *r = (float*) malloc(sizeof(float) * width * height);
    float *g = (float*) malloc(sizeof(float) * width * height);
    float *b = (float*) malloc(sizeof(float) * width * height);
    return new_ImageSoA1(width, height, r, g, b);
}

void image_delete(ImageSoA_t* img) {
    if (img != NULL) {
        if (image_getR(img) != NULL) {
            free(image_getR(img));
        }
        if (image_getG(img) != NULL) {
            free(image_getG(img));
        }
        if (image_getB(img) != NULL) {
            free(image_getB(img));
        }
        free(img);
    }
}


void image_setPixel(ImageSoA_t* img, int x, int y, int c, float val) {
    int width = image_getWidth(img);

    float *r = image_getR(img);
    float *g = image_getG(img);
    float *b = image_getB(img);

    switch (c) {
        case 0: {
            r[y * width + x] = val;
            break;
        }
        case 1: {
            g[y * width + x] = val;
            break;
        }
        case 2: {
            b[y * width + x] = val;
            break;
        }
        default:
            std::cerr << "Wrong channel, it must be 0, 1 or 2" << std::endl;
    }

    return;
}

float image_getPixel(ImageSoA_t* img, int x, int y, int c) {
    int width = image_getWidth(img);

    float *r = image_getR(img);
    float *g = image_getG(img);
    float *b = image_getB(img);

    switch (c) {
        case 0: {
            return r[y * width + x];
        }
        case 1: {
            return g[y * width + x];
        }
        case 2: {
            return b[y * width + x];
        }
        default:
            std::cerr << "Wrong channel, it must be 0, 1 or 2" << std::endl;
    }
    return 0;
}