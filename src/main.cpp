#include "image/PpmParser.h"
#include "kernel/Kernel.h"
#include "utils/Utils.h"

int main() {
//    test_images();

    Image_t* image = PPM_import("../resources/computer_programming.ppm");
    int width = image->width;
    int height = image->height;
    int channels = image->channels;

    float* kernel = createKernel(kernelsType::gaussianBlur);
    int kernelRadius = kernelWidth() / 2;
    int kernelW = kernelWidth();

    float* imageIter = image->data;
    float* kernelIter = kernel;

    Image_t* result = new_image(width, height, channels);
    float* resultIter = result->data;

    float imagePixel;
    float kernelValue;
    float accum;
    int xOffset;
    int yOffset;

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
                            kernelValue = kernelIter[(y + kernelRadius) * kernelW + x + kernelRadius];
                            accum += (imagePixel * kernelValue);
                        }
                    }
                }
            resultIter[(i * width + j) * channels + k] = clamp(accum, 0, 1);
            }
        }
    }

    PPM_export("../resources/results/output.ppm", result);

    free(kernel);
    image_delete(image);
    image_delete(result);
    return 0;
}

