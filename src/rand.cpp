#include "rand.h"

#include <cstdlib>

void rand_seed(int seed) {
    std::srand(static_cast<unsigned int>(seed));
}

