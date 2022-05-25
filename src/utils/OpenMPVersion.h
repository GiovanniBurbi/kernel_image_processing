//
// Created by giovanni on 16/05/22.
//

#ifndef KERNEL_IMAGE_PROCESSING_OPENMPVERSION_H
#define KERNEL_IMAGE_PROCESSING_OPENMPVERSION_H
#include <unordered_map>
#include <iostream>
#include <omp.h>

// simd directive is supported from OpenMP 4.0
static inline void checkVersionOpenMP(){
    std::unordered_map<unsigned,std::string> map{
            {200505,"2.5"},{200805,"3.0"},{201107,"3.1"},{201307,"4.0"},{201511,"4.5"},{201811,"5.0"},{202011,"5.1"}};
    std::cout << "We have OpenMP " << map.at(_OPENMP) << ".\n";
}
#endif //KERNEL_IMAGE_PROCESSING_OPENMPVERSION_H
