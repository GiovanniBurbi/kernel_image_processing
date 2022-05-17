#include "image/PpmParser.h"
#include "kernel/Kernel.h"
#include "convolution/Convolution.h"
#include "utils/OpenMPVersion.h"

#include <chrono>

#define IMPORT_PATH "../resources/"
#define EXPORT_PATH "../resources/results/"
#define IMAGE "tiger"

#define SEQUENTIAL true
#define PARALLEL false
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
        std::string log;
        log.append("Sequential version ");

        std::chrono::high_resolution_clock::time_point startTime;
        std::chrono::high_resolution_clock::time_point endTime;

        if (UNROLLING) {
            log.append("with unrolling ");
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
            log.append("naive ");
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
        log.append("took ").append(std::to_string(time/ITER)).append(" seconds");
        printf("%s\n", log.c_str());
    }

    if (PARALLEL) {
        std::string log;
        log.append("OpenMP parallel version ");

        std::chrono::high_resolution_clock::time_point startTime;
        std::chrono::high_resolution_clock::time_point endTime;

        int nThreads = 2;

        if(UNROLLING) {
            log.append("with unrolling ");

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
            log.append("naive ");

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
        log.append("took ").append(std::to_string(time/ITER)).append(" seconds");
        printf("%s\n", log.c_str());
    }

    outputname.append(".ppm");
    PPM_export(outputname.c_str(), output);

    free(kernel);
    image_delete(image);
    image_delete(output);

    return 0;
}

