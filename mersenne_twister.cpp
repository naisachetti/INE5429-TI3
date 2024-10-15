// Implementação baseada na descrição do algoritmo de https://www.sciencedirect.com/topics/computer-science/mersenne-twister.

#include <iostream>

#define N 624
#define M 397
#define UPPER_MASK 0x80000000 // w-r bits mais significativos
#define LOWER_MASK 0x7fffffff // r bits menos significativos
static unsigned long state[N]; // Vetor de estados Y
static int stateIdx = N+1; // stateIdx==N+1 significa que o vetor de estados não foi inicializado

uint32_t genrand_int32(void) {
    static uint32_t constants[2] = {0x0, 0x9908b0df};
    uint32_t y;
    // Atualização do vetor de estados
    if (stateIdx >= N) {
        int k;
        for (k=0; k<N-M; k++) {
            y = (state[k]&UPPER_MASK) | (state[k+1]&LOWER_MASK);
            state[k] = state[k+M] ^ (y>>1) ^ constants[y&0x1];
        }
        for (; k<N-1; k++) {
            y = (state[k]&UPPER_MASK) | (state[k+1]&LOWER_MASK);
            state[k] = state[k+(M-N)] ^ (y>>1) ^ constants[y&0x1];
        }
        y = (state[N-1]&UPPER_MASK) | (state[0]&LOWER_MASK);
        state[N-1] = state[M-1] ^ (y>>1) ^ constants[y&0x1];
        stateIdx = 0;
    }
    // Geração de número pseudoaleatório
    y = state[stateIdx++];
    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= (y >> 18);

    return y;
}

int main() {
    // Inicialização do vetor de estados
    for (int k = 0; k < N; k++) {
        state[k] = k;
    }

    // Geração de números pseudoaleatórios
    for (int k = 0; k < 10000; k++) {
        printf("%u\n", genrand_int32());
    }
}