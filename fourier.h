#ifndef FOURIER_H
#define FOURIER_H


extern char hdm8[8][8],hdm16[16][16], hdm32[32][32], hdm64[64][64], hdm128[128][128], hdm256[256][256];

extern void initialize_fourier();

extern double fourier_coef[65535];
extern double influence[15];
extern double coupling[15];
extern double total_influence;

extern void calculate_fourier(unsigned int x,bool (*bool_func)(unsigned int,unsigned int));
extern double spectral_entropy();
extern double calculate_influence(unsigned int j);
extern double calculate_influence2(unsigned int j);
extern double calculate_coupling(unsigned int i, unsigned int j,bool (*bool_func)(unsigned int,unsigned int));
extern double calculate_coupling2(unsigned int i, unsigned int j);
extern double total_coupling(unsigned int i);
extern double max(double a, double b);
extern unsigned int degree();
extern double levelWeight(unsigned k);
extern double Stab(double r);
extern double NoiseSensitivity(double r);

extern double C;

#endif // FOURIER_H
