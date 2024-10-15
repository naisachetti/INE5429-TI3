#include "mersenne_twister_big.h"
#include <chrono>
#include <iostream>
#include <gmp.h>

int main() {
    MyMersenneTwister mt(4096,1); // MyMersenneTwister(uint32_t num_bits, uint32_t seed);

    // Geração de números pseudoaleatórios
    mpz_t prn;
    mpz_init(prn);
    double total_time = 0;
    for (uint32_t k = 0; k < 1000000; k++) {
        auto start = std::chrono::steady_clock::now();
        mt.genrand_int(prn);
        auto end = std::chrono::steady_clock::now();
        total_time += std::chrono::duration<double, std::micro>(end-start).count();

        std::cout << prn << std::endl;
    }
    // std::cout << "Tempo médio de geração de cada número: " << (total_time/1000000) << " µs" << std::endl;
    mpz_clear(prn);
}