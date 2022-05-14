#include "image/PpmParser.h"

#include <iostream>

int main() {
//    test_images();

    Image_t* inputImg = PPM_import("../resources/computer_programming.ppm");

    return 0;
}
