// Implementação baseada na descrição do algoritmo de https://www.sciencedirect.com/topics/computer-science/mersenne-twister.

#ifndef __MERSENNE_TWISTER_BIG_H
#define __MERSENNE_TWISTER_BIG_H

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <gmp.h>

#define N 624
#define M 397

class MyMersenneTwister {
    public:
    MyMersenneTwister(uint32_t num_bits, uint32_t seed);
    ~MyMersenneTwister();
    void genrand_int(mpz_t &y);

    private:
    uint32_t num_bits;
    uint32_t seed;
    mpz_t state[N]; // Vetor de estados
    int stateIdx = N+1; // stateIdx==N+1 significa que o vetor de estados não foi inicializados
    // Variáveis sempre utilizadas pela função genrand_int()
    mpz_t upper_mask, lower_mask;
    mpz_t constants[2];
    mpz_t aux1, aux2, aux3;
    mpz_t shift1, shift2;
    mpz_t bitmask_numbits;
};

#endif /* #ifndef __MTWISTER_H */