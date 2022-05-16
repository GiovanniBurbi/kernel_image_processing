#include "image/PpmParser.h"
#include "kernel/Kernel.h"
#include "convolution/Convolution.h"
#include "utils/OpenMPVersion.h"

#include <omp.h>
#include <chrono>

#define SEQUENTIAL true
#define PARALLEL false
#define UNROLLING false


int main() {
//    checkVersionOpenMP();

//    Image_t* image = PPM_importWithPadding("../resources/tiger.ppm");

    Image_t* image = PPM_import("../resources/tiger.ppm");
    float* kernel = createKernel(kernelsType::outline);

    Image_t* result;

    if (SEQUENTIAL) {
        std::chrono::high_resolution_clock::time_point startTime;
        std::chrono::high_resolution_clock::time_point endTime;
        if (UNROLLING) {
            startTime = std::chrono::high_resolution_clock::now();
            result = convolutionUnrolling(image, kernel, kernelWidth());
            endTime = std::chrono::high_resolution_clock::now();
        } else {
            startTime = std::chrono::high_resolution_clock::now();
            result = convolution(image, kernel, kernelWidth());
            endTime = std::chrono::high_resolution_clock::now();

        }
        float time = std::chrono::duration_cast<std::chrono::duration<float>>(endTime - startTime).count();
        printf("sequential time %f", time);
    }

    if (PARALLEL) {
        std::chrono::high_resolution_clock::time_point startTime;
        std::chrono::high_resolution_clock::time_point endTime;

        int nThreads = 2;

        if(UNROLLING) {
            startTime = std::chrono::high_resolution_clock::now();
            result = convolutionOMPUnrollingDoubleSIMD(image, kernel, kernelWidth(), nThreads);
            endTime = std::chrono::high_resolution_clock::now();
        } else {
            startTime = std::chrono::high_resolution_clock::now();
            result = convolutionOMPNaiveNoBorder(image, kernel, kernelWidth(), nThreads);
            endTime = std::chrono::high_resolution_clock::now();
        }
        float time = std::chrono::duration_cast<std::chrono::duration<float>>(endTime - startTime).count();
        printf("sequential time %f", time);
    }

    PPM_export("../resources/results/output8.ppm", result);

    free(kernel);
    image_delete(image);
    image_delete(result);
    return 0;
}

