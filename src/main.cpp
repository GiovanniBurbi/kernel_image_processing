#include "image/PpmParser.h"
#include "kernel/Kernel.h"
#include "convolution/Convolution.h"
#include "utils/OpenMPVersion.h"

#include <chrono>

#define IMPORT_PATH "../resources/"
#define EXPORT_PATH "../resources/results/"
#define IMAGE "tiger"

#define SEQUENTIAL false
#define PARALLEL true
#define UNROLLING false
#define ITER 1


int main() {
//    checkVersionOpenMP();

    std::string filename;
    std::string outputname;

    filename.append(IMPORT_PATH).append(IMAGE).append(".ppm");
    outputname.append(EXPORT_PATH).append(IMAGE);

//    Image_t* image = PPM_importWithPadding("../resources/tiger.ppm");

    Image_t* image = PPM_import(filename.c_str());
    float* kernel = createKernel(kernelsType::emboss);

    Image_t* output;
    float time = 0;

    if (SEQUENTIAL) {
        std::chrono::high_resolution_clock::time_point startTime;
        std::chrono::high_resolution_clock::time_point endTime;

        if (UNROLLING) {
            for (int i = 0; i < ITER; i++) {
                startTime = std::chrono::high_resolution_clock::now();
                Image_t *result = convolutionUnrolling(image, kernel);
                endTime = std::chrono::high_resolution_clock::now();
                time += std::chrono::duration_cast<std::chrono::duration<float>>(endTime - startTime).count();
                if (i != ITER - 1)
                    image_delete(result);
                else output = result;
            }
            outputname.append("Sequential").append("Unrolling");
        } else {
            for (int i = 0; i < ITER; i++) {
                startTime = std::chrono::high_resolution_clock::now();
                Image_t *result = convolution(image, kernel);
                endTime = std::chrono::high_resolution_clock::now();
                time += std::chrono::duration_cast<std::chrono::duration<float>>(endTime - startTime).count();
                if (i != ITER - 1)
                    image_delete(result);
                else output = result;
            }
            outputname.append("Sequential").append("Naive");
        }

        float time = std::chrono::duration_cast<std::chrono::duration<float>>(endTime - startTime).count();
        printf("Sequential time %f", time / ITER);
    }

    if (PARALLEL) {
        std::chrono::high_resolution_clock::time_point startTime;
        std::chrono::high_resolution_clock::time_point endTime;

        int nThreads = 2;

        if(UNROLLING) {
            for (int i = 0; i < ITER; i++) {
                startTime = std::chrono::high_resolution_clock::now();
//            Image_t* result = convolutionOMPUnrolling(image, kernel, nThreads);
                Image_t* result = convolutionOMPUnrollingDoubleSIMD(image, kernel, nThreads);
                endTime = std::chrono::high_resolution_clock::now();
                time += std::chrono::duration_cast<std::chrono::duration<float>>(endTime - startTime).count();
                if (i != ITER - 1)
                    image_delete(result);
                else output = result;
            }
            outputname.append("Openmp").append("Parallel").append("Unrolling").append(std::to_string(nThreads)).append("Threads");
        } else {
            for (int i = 0; i < ITER; i++) {
                startTime = std::chrono::high_resolution_clock::now();
                Image_t *result = convolutionOMPNaive(image, kernel, nThreads);
                endTime = std::chrono::high_resolution_clock::now();
                time += std::chrono::duration_cast<std::chrono::duration<float>>(endTime - startTime).count();
                if (i != ITER - 1)
                    image_delete(result);
                else output = result;
            }
            outputname.append("Openmp").append("Parallel").append("Naive").append(std::to_string(nThreads)).append("Threads");
        }

        printf("Parallel time %f", time / ITER);
    }

    outputname.append(".ppm");
    PPM_export(outputname.c_str(), output);

    free(kernel);
    image_delete(image);
    image_delete(output);

    return 0;
}

