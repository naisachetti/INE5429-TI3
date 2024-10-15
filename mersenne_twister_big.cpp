// Implementação baseada na descrição do algoritmo de https://www.sciencedirect.com/topics/computer-science/mersenne-twister.

#include "mersenne_twister_big.h"

MyMersenneTwister::MyMersenneTwister(uint32_t num_bits, uint32_t seed) {
    this->num_bits = num_bits;
    this->seed = seed;

    /*
    * Definição dos números hexadecimais usados na geração de PRNs.
    *
    * Optou-se por sempre utilizar números com incrementos de 32 bits para que fosse possível replicar as constantes do algoritmo original, 
    * de acordo com a necessidade. Por exemplo: Para a geração de números de 32 bits, utiliza-se a constante "0x9d2c5680" para relizar o 
    * primeiro shift da fase de tempering do Mersenne Twister. Já para 33 bits até 64 bits, utiliza-se "0x9d2c5680 9d2c5680", e assim por diante.
    * 
    * Também optou-se por iniciar os valores do vetor de estados (state) em "0x1000 0000" para núeros de 1 a 32 bits, "0x1000 0000 0000 0000" 
    * para números de 33 a 64 bits, e assim por diante (ver init_state_num), para possibilitar números com grandes quantidades de bits fossem 
    * gerados. Os valores de state são definidos como incrementos de 1 sucessivos sobre "0x1000 0000" e embaralhamentos feitos com o uso de uma semente
    * (ver inicialização do vetor de estados, nesta função).
    */
    std::string upper_mask_hex = "8000 0000"; // Usado ao inicializar upper_mask. Na descrição original do algoritmo, os w-r bits mais significativos, utilizando w-r = 1.
    std::string lower_mask_hex = "7fff ffff"; // Usado ao inicializar upper_mask. Na descrição original do algoritmo, os r bits menos significativos.
    std::string constants_i1 = "9908 b0df"; // Usado ao inicializar constants[1]
    std::string shift1_hex = "9d2c 5680"; // Usado ao inicializar shift1
    std::string shift2_hex = "efc6 0000"; // Usado ao inicializar shift2
    std::string init_state_num = "1000 0000"; // Usado ao inicializar state.
    for (int k=0; k<std::ceil(float(num_bits)/32)-1; k++) {
        upper_mask_hex += " 0000 0000";
        lower_mask_hex += " ffff ffff";
        constants_i1 += " 9908 b0df";
        shift1_hex += " 9d2c 5680";
        shift2_hex += " efc6 0000";
        init_state_num += "0000 0000";
    }

    mpz_init(this->upper_mask);
    mpz_init(this->lower_mask);
    mpz_set_str(this->upper_mask, &upper_mask_hex[0], 16);
    mpz_set_str(this->lower_mask, &lower_mask_hex[0], 16);

    mpz_init(this->constants[0]);
    mpz_init(this->constants[1]);
    mpz_set_ui(this->constants[0], 0);
    mpz_set_str(this->constants[1], &constants_i1[0], 16);

    mpz_init(this->aux1);
    mpz_init(this->aux2);
    mpz_init(this->aux3);

    mpz_init(this->shift1);
    mpz_init(this->shift2);
    mpz_set_str(this->shift1, &shift1_hex[0], 16);
    mpz_set_str(this->shift2, &shift2_hex[0], 16);

    // Inicialização do vetor de estados.
    mpz_t state_num;
    mpz_init(state_num);
    mpz_set_str(state_num, &init_state_num[0], 16);
    mpz_t seed_particle;
    mpz_init(seed_particle);
    mpz_set_ui(seed_particle, this->seed);
    uint32_t porpoc_num_bits = std::ceil(num_bits/32); // Utilizado para preparar seed_particle proporcionalmente ao número de bits
    for (int k = 0; k < N; k++) {
        // Alterando o estado com state_num, para permitir a geração de números grandes
        mpz_init(this->state[k]);
        mpz_set(this->state[k], state_num); // state[k] = init_state_num + k
        mpz_add_ui(state_num, state_num, 1); // state_num += 1

        // Preparando seed_particle (operações escolhidas de modo a tentar diferenciar dígitos de mais alto valor dos números entre 
        // uma semente e outra, de acordo com o tamanho do número).
        mpz_set(this->aux1, seed_particle); // aux1 = seed_particle
        mpz_mul_2exp(seed_particle, seed_particle, 12*porpoc_num_bits); // seed_particle <<= 12*porpoc_num_bits
        mpz_add_ui(seed_particle, seed_particle, 3*porpoc_num_bits); // seed_particle += 3*porpoc_num_bits
        mpz_mul_ui(seed_particle, seed_particle, 7*porpoc_num_bits); // seed_particle *= 7*porpoc_num_bits
        mpz_mul_2exp(seed_particle, seed_particle, 5*porpoc_num_bits); // seed_particle <<= 5*porpoc_num_bits
        mpz_xor(seed_particle, seed_particle, this->aux1); // seed_particle ^= aux1

        // Alterando o estado com seed_particle, para que o gerador gere números com base em uma semente
        mpz_xor(state[k], state[k], seed_particle); // state[k] ^= seed_particle
    }
    mpz_clear(state_num);

    // Inicialização da máscara para a geração de números com o número de bits escolhido
    std::string bitmask_numbits_str = "";
    for (int k=0; k<num_bits; k++) {
        bitmask_numbits_str += "1";
    }
    mpz_init(this->bitmask_numbits);
    mpz_set_str(this->bitmask_numbits, &bitmask_numbits_str[0], 2);
}

MyMersenneTwister::~MyMersenneTwister() {
    mpz_clear(upper_mask);
    mpz_clear(lower_mask);
    mpz_clear(constants[0]);
    mpz_clear(constants[1]);
    mpz_clear(aux1);
    mpz_clear(aux2);
    mpz_clear(aux3);
    mpz_clear(shift1);
    mpz_clear(shift2);
    mpz_clear(bitmask_numbits);
    for (int k=0; k < N; k++) {
        mpz_clear(state[k]);
    }
}

void MyMersenneTwister::genrand_int(mpz_t &y) {
    // Atualização do vetor de estados
    if (this->stateIdx >= N) {
        int k;
        for (k=0; k<N-M; k++) {
            // Equivalente a: y = (state[k]&UPPER_MASK) | (state[k+1]&LOWER_MASK);
            mpz_and(this->aux1, this->state[k], this->upper_mask);
            mpz_and(this->aux2, this->state[k+1], this->lower_mask);
            mpz_ior(y, this->aux1, this->aux2);
            // state[k] = state[k+M] ^ (y>>1) ^ constants[y&0x1];
            mpz_tdiv_q_2exp(this->aux1, y, 1); // y>>1
            mpz_xor(this->aux2, this->state[k+M], this->aux1); // state[k+M] ^ (y>>1)
            mpz_set_ui(this->aux1, 1); 
            mpz_and(this->aux3, y, this->aux1); // y&0x1
            mpz_xor(this->state[k], this->aux2, this->constants[mpz_get_ui(this->aux3)]); // state[k] = state[k+M] ^ (y>>1) ^ constants[y&0x1];
        }
        for (; k<N-1; k++) {
            // y = (state[k]&UPPER_MASK) | (state[k+1]&LOWER_MASK);
            mpz_and(this->aux1, this->state[k], this->upper_mask);
            mpz_and(this->aux2, this->state[k+1], this->lower_mask);
            mpz_ior(y, this->aux1, this->aux2);
            // state[k] = state[k+(M-N)] ^ (y>>1) ^ constants[y&0x1];
            mpz_tdiv_q_2exp(aux1, y, 1); // y>>1
            mpz_xor(this->aux2, this->state[k+(M-N)], this->aux1); // state[k+(M-N)] ^ (y>>1)
            mpz_set_ui(this->aux1, 1); 
            mpz_and(this->aux3, y, aux1); // y&0x1
            mpz_xor(this->state[k], this->aux2, this->constants[mpz_get_ui(this->aux3)]); // state[k] = state[k+(M-N)] ^ (y>>1) ^ constants[y&0x1];
        }
        // y = (state[N-1]&UPPER_MASK) | (state[0]&LOWER_MASK);
        mpz_and(this->aux1, this->state[N-1], this->upper_mask);
        mpz_and(this->aux2, this->state[0], this->lower_mask);
        mpz_ior(y, this->aux1, this->aux2);
        // state[N-1] = state[M-1] ^ (y>>1) ^ constants[y&0x1];
        mpz_tdiv_q_2exp(this->aux1, y, 1); // y>>1
        mpz_xor(this->aux2, this->state[M-1], this->aux1); // state[M-1] ^ (y>>1)
        mpz_set_ui(this->aux1, 1); 
        mpz_and(this->aux3, y, this->aux1); // y&0x1
        mpz_xor(this->state[N-1], this->aux2, this->constants[mpz_get_ui(this->aux3)]); // state[N-1] = state[M-1] ^ (y>>1) ^ constants[y&0x1];

        this->stateIdx = 0;
    }
    // Geração de número pseudoaleatório
    // y = state[stateIdx++];
    mpz_set(y, this->state[this->stateIdx++]);
    // Tempering
    // y ^= (y >> 11);
    mpz_tdiv_q_2exp(this->aux1, y, 11);
    mpz_xor(y, y, this->aux1);
    // y ^= (y << 7) & 0x9d2c5680;
    mpz_mul_2exp(this->aux1, y, 7);
    mpz_and(this->aux3, this->aux1, this->shift1);
    mpz_xor(y, y, this->aux3);
    // y ^= (y << 15) & 0xefc60000;
    mpz_mul_2exp(this->aux1, y, 15);
    mpz_and(this->aux3, this->aux1, this->shift2);
    mpz_xor(y, y, this->aux3);

    // Aplicação da máscara para a geração de números com o número de bits escolhido
    mpz_and(y, y, this->bitmask_numbits);

    // y ^= (y >> 18);
    mpz_tdiv_q_2exp(this->aux1, y, 18);
    mpz_xor(y, y, this->aux1);  
}