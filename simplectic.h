#include <stdint.h>

#define TABLE_SIZE 256

typedef struct simplectic_seed {
    uint8_t table[TABLE_SIZE * 2];
} simplectic_seed;

/**
 * Create a new seed. Will write new pseudo-random values to seed, for
 * possible re-use for creating new seeds in the future.
 */
void simplectic_seed_fill(simplectic_seed *sd, uint32_t seed);

double simplectic2(const double point[const 2], const simplectic_seed *seed);
double simplectic3(const double point[const 3], const simplectic_seed *seed);
double simplectic4(const double point[const 4], const simplectic_seed *seed);
//double simplecticn(const simplectic_pointn *point, ssize_t num_dims, const simplectic_seed *seed);
