#include <cassert>
#include <chrono>

#include "image/PpmParser.h"
#include "kernel/Kernel.h"
#include "convolution/Convolution.h"
#include "utils/OpenMPVersion.h"

#define IMPORT_PATH "../resources/source/"
#define EXPORT_PATH "../resources/results/"
#define IMAGE "lake"

#define SEQUENTIAL false
#define PARALLEL true
#define UNROLLING false
#define SOA false
#define ITER 15
#define N_THREADS 12


int main() {
//    checkVersionOpenMP();
    assert((SEQUENTIAL && PARALLEL) != true);

    std::string filename;
    std::string output_name;

    filename.append(IMPORT_PATH).append(IMAGE).append(".ppm");
    output_name.append(EXPORT_PATH).append(IMAGE);

    auto maskType = kernelsType::outline;

    float* kernel = createKernel(maskType);
    std::string maskName = kernelName(maskType);

    if (!SOA) {
        float time = 0;

        Image_t* image = PPM_import(filename.c_str());
        Image_t* output;

        if (SEQUENTIAL) {
            std::string log;
            log.append("Sequential version ");

            std::chrono::high_resolution_clock::time_point startTime;
            std::chrono::high_resolution_clock::time_point endTime;

            if (UNROLLING) {
                log.append("with unrolling ");
                for (int i = 0; i < ITER; i++) {
                    startTime = std::chrono::high_resolution_clock::now();
//                    For RGB image only
//                    Image_t *result = convolutionUnrollingChannels(image, kernel);
//                    For RGB image and 3x3 kernels
                    Image_t *result = convolutionUnrolling(image, kernel);
                    endTime = std::chrono::high_resolution_clock::now();
                    time += std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                    if (i != ITER - 1)
                        image_delete(result);
                    else output = result;
                }
                output_name.append(maskName).append("Sequential").append("Unrolling");
            } else {
                log.append("naive ");
                for (int i = 0; i < ITER; i++) {
                    startTime = std::chrono::high_resolution_clock::now();
//                    Naive generic convolution
                    Image_t *result = convolution(image, kernel);
                    endTime = std::chrono::high_resolution_clock::now();
                    time += std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                    if (i != ITER - 1)
                        image_delete(result);
                    else output = result;
                }
                output_name.append(maskName).append("Sequential").append("Naive");
            }

            log.append("took ").append(std::to_string(time/ITER)).append(" microseconds");
            printf("%s\n", log.c_str());
        }

        if (PARALLEL) {
            std::string log;
            log.append("OpenMP parallel version ").append("with ").append(std::to_string(N_THREADS)).append(" threads ");

            std::chrono::high_resolution_clock::time_point startTime;
            std::chrono::high_resolution_clock::time_point endTime;

            if(UNROLLING) {
                log.append("with unrolling ");

                for (int i = 0; i < ITER; i++) {
                    startTime = std::chrono::high_resolution_clock::now();
//                   For 3x3 kernels only
                    Image_t* result = convolutionOMPUnrollingSIMDChannels(image, kernel, N_THREADS);
//                  For RGB image and 3x3 kernel only
//                    Image_t* result = convolutionOMPUnrolling(image, kernel, N_THREADS);
//                  For RGB image only
//                    Image_t *result = convolutionOMPUnrollingChannels(image, kernel, N_THREADS);
                    endTime = std::chrono::high_resolution_clock::now();
                    time += std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                    if (i != ITER - 1)
                        image_delete(result);
                    else output = result;
                }
                output_name.append(maskName).append("Openmp").append("Parallel").append("Unrolling").append(std::to_string(N_THREADS)).append("Threads");
            } else {
                log.append("naive ");

                for (int i = 0; i < ITER; i++) {
                    startTime = std::chrono::high_resolution_clock::now();
//                    Generic naive convolution
                    Image_t *result = convolutionOMPNaive(image, kernel, N_THREADS);
                    endTime = std::chrono::high_resolution_clock::now();
                    time += std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                    if (i != ITER - 1)
                        image_delete(result);
                    else output = result;
                }
                output_name.append(maskName).append("Openmp").append("Parallel").append("Naive").append(std::to_string(N_THREADS)).append("Threads");
            }

            log.append("took ").append(std::to_string(time/ITER)).append(" microseconds");
            printf("%s\n", log.c_str());
        }

        output_name.append(".ppm");

        PPM_export(output_name.c_str(), output);

        image_delete(image);
        image_delete(output);
    }

    if (SOA) {
        float time = 0;

        ImageSoA_t* image = PPM_importSoA(filename.c_str());
        ImageSoA_t* output;

        if (SEQUENTIAL) {
            std::string log;
            log.append("Sequential version ");

            std::chrono::high_resolution_clock::time_point startTime;
            std::chrono::high_resolution_clock::time_point endTime;

            if (UNROLLING) {
                log.append("with unrolling ");
                for (int i = 0; i < ITER; i++) {
                    startTime = std::chrono::high_resolution_clock::now();
//                    For RGB image and 3x3 kernels only
                    ImageSoA_t *result = convolutionUnrollingSoA(image, kernel);
                    endTime = std::chrono::high_resolution_clock::now();
                    time += std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                    if (i != ITER - 1)
                        image_delete(result);
                    else output = result;
                }
                output_name.append(maskName).append("Sequential").append("Unrolling");
            } else {
                log.append("naive ");
                for (int i = 0; i < ITER; i++) {
                    startTime = std::chrono::high_resolution_clock::now();
//                    For RGB image only
                    ImageSoA_t *result = convolutionSoA(image, kernel);
                    endTime = std::chrono::high_resolution_clock::now();
                    time += std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                    if (i != ITER - 1)
                        image_delete(result);
                    else output = result;
                }
                output_name.append(maskName).append("Sequential").append("SoA").append("Naive");
            }

            log.append("took ").append(std::to_string(time/ITER)).append(" microseconds");
            printf("%s\n", log.c_str());
        }

        if (PARALLEL) {
            std::string log;
            log.append("OpenMP parallel version ").append("with ").append(std::to_string(N_THREADS)).append(" threads ");

            std::chrono::high_resolution_clock::time_point startTime;
            std::chrono::high_resolution_clock::time_point endTime;

            if(UNROLLING) {
                log.append("with unrolling ");

                for (int i = 0; i < ITER; i++) {
                    startTime = std::chrono::high_resolution_clock::now();
//                    For RGB image and 3x3 kernels only
                    ImageSoA_t* result = convolutionOMPUnrollingSoA(image, kernel, N_THREADS);
                    endTime = std::chrono::high_resolution_clock::now();
                    time += std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                    if (i != ITER - 1)
                        image_delete(result);
                    else output = result;
                }
                output_name.append(maskName).append("Openmp").append("Parallel").append("Unrolling").append(std::to_string(N_THREADS)).append("Threads");
            } else {
                log.append("naive ");

                for (int i = 0; i < ITER; i++) {
                    startTime = std::chrono::high_resolution_clock::now();
//                    For RGB images only
                    ImageSoA_t *result = convolutionOMPNaiveSoA(image, kernel, N_THREADS);
                    endTime = std::chrono::high_resolution_clock::now();
                    time += std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
                    if (i != ITER - 1)
                        image_delete(result);
                    else output = result;
                }
                output_name.append(maskName).append("Openmp").append("Parallel").append("Naive").append(std::to_string(N_THREADS)).append("Threads");
            }

            log.append("took ").append(std::to_string(time/ITER)).append(" microseconds");
            printf("%s\n", log.c_str());
        }

        output_name.append(".ppm");

        PPM_exportSoA(output_name.c_str(), output);

        image_delete(image);
        image_delete(output);
    }

    free(kernel);

    return 0;
}