cmake_minimum_required(VERSION 3.21)
project(kernel_image_processing)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(kernel_image_processing src/main.cpp src/image/Image.cpp src/image/Image.h src/utils/Utils.h src/image/PpmParser.cpp src/image/PpmParser.h src/kernel/Kernel.cpp src/kernel/Kernel.h src/convolution/Convolution.cpp src/convolution/Convolution.h src/utils/OpenMPVersion.h src/image/ImageSoA.cpp src/image/ImageSoA.h)
include_directories(./src)
