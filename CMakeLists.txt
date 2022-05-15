cmake_minimum_required(VERSION 3.21)
project(kernel_image_processing)

set(CMAKE_CXX_STANDARD 14)

add_executable(kernel_image_processing src/main.cpp src/image/Image.cpp src/image/Image.h src/utils/Utils.h src/image/PpmParser.cpp src/image/PpmParser.h src/kernel/types/BoxBlur.h src/kernel/types/GaussianBlur.h src/kernel/types/Emboss.h src/kernel/types/Outline.h src/kernel/types/Sharpen.h src/kernel/Kernel.cpp src/kernel/Kernel.h)
include_directories(./src)
