//
// Created by giovanni on 13/05/22.
//

#ifndef KERNEL_IMAGE_PROCESSING_PPMPARSER_H
#define KERNEL_IMAGE_PROCESSING_PPMPARSER_H

#include "Image.h"
#include "ImageSoA.h"

Image_t* PPM_import(const char *filename);
ImageSoA_t* PPM_importSoA(const char *filename);
Image_t* PPM_importWithPadding(const char *filename);
bool PPM_export(const char *filename, Image_t* img);
bool PPM_exportSoA(const char *filename, ImageSoA_t* img);

void test_images();

#endif //KERNEL_IMAGE_PROCESSING_PPMPARSER_H
