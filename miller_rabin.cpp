// Código implementado de acordo com a descrição do algoritmo explicada em aula

#include <iostream>
#include <stdio.h>
#include <bitset>
#include <random>
#include <cmath>
#include <float.h>

uint32_t modular_exp(uint32_t base, uint32_t exp, uint32_t mod_num) {
    // Faz a exponenciação modular de acordo com o algoritmo descrito em https://people.reed.edu/~jerry/361/lectures/bigprimes.pdf
    uint32_t x = 1;
    uint32_t y = base;
    uint32_t f = exp;
    while (f > 0) {
        if (f % 2 == 0) {
            y = (y*y) % mod_num;
            f /= 2;
        } else {
            x = (x*y) % mod_num;
            f -= 1;
        }
    }
    return x;
}

bool miller_rabin(uint32_t num, uint32_t num_it) {
    // Escreva n-1=(2^k)*m, onde m é ímpar
    // Obtendo k e m via contagem de bits 0 menos significativos em n-1
    std::bitset<32> bs(num-1);
    uint32_t k = 0;
    uint32_t m = num-1;
    for (int i = 0; i < bs.size(); i++) {
        if (bs[i] == 0) {
            k++;
            m >>= 1;
        } else break;
    }

    std::random_device device;
    std::mt19937 mt(device()); 
    std::uniform_int_distribution<std::mt19937::result_type> dist(2,num-1);
    // Repete para um número num_it de a's aleatórios, sendo 1< a < n.
    for (uint32_t it=0; it<num_it; it++) {
        uint32_t a = dist(mt);
        // Se (a^m) mod n = 1
        if (modular_exp(a,m,num) == 1) {
            return true;
        }

        // Para i=0 até k-1
        for (int i=0; i<=k-1; i++) {
            // Se a^((2^i)*m) mod n = n-1
            uint32_t exponent = std::pow(2,i)*m;
            if (modular_exp(a,exponent,num) == num-1) {
                return true;
            }
        }
    }

    return false;
}

// Exemplo de chamada: ./millerrabin num_it
int main(int argc, char *argv[]) {
    uint32_t num_it = std::atol(argv[2]); // Quantidade de iterações para testar se um número é primo (isso aumenta a probabilidade 
                                              // do número realmente ser primo, se passar no teste).

    std::random_device device;
    std::mt19937 mt(device()); 
    std::uniform_int_distribution<std::mt19937::result_type> dist(2,UINT32_MAX);
    while (true) {
        uint32_t num_test = dist(mt);
        if (num_test % 2 == 0) {
            num_test -= 1;
        }
        bool passed = miller_rabin(num_test, num_it);
        if (passed) {
            printf("%u\n", num_test);
            return 0;
        }
    }
}