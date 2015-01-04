#include "simplectic.h"
#include <stdio.h>
#include <math.h>

#define CL(a, b, c) (a < b ? b : a > c ? c : a)

int main(int argc, char **argv) {
    // generate a 1024x1024 PGM image
    simplectic_seed sd;
    simplectic_seed_fill(&sd, 0);
    puts("P2\n1024 1024\n255");
    double point[2];
    double val;

    for (int y = 0; y < 256; y++) {
        for (int x = 0; x < 256; x++) {
            point[0] = x / 64.0;
            point[1] = y / 64.0;
            val = simplectic2(point, &sd);
            val = CL(val * 0.5 + 0.5, 0.0, 1.0) * 255.0;
            printf("%d ", (uint8_t)val);
        }
        puts("");
    }
    return 0;
}
