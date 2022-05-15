#include "image/PpmParser.h"
#include "kernel/Kernel.h"
#include "convolution/Convolution.h"

int main() {
    Image_t* image = PPM_import("../resources/computer_programming.ppm");

    float* kernel = createKernel(kernelsType::emboss);

    Image_t* result = new_image(image->width, image->height, image->channels);

    convolution(image, result, kernel, kernelWidth());

    PPM_export("../resources/results/output.ppm", result);

    free(kernel);
    image_delete(image);
    image_delete(result);
    return 0;
}

