#ifndef BOOLEAN_H
#define BOOLEAN_H
#include<stdio.h>
#define AND 0
#define OR  1
#define INPUTX 8
extern unsigned INPUT; //MAX 8


extern unsigned fact(unsigned int n);
extern void print_binary(unsigned int k);
extern void fprint_binary(FILE *fp, unsigned int k);
extern float correlation(float x[], float y[]);
extern bool boolean_generic(unsigned int no, unsigned int input);
extern bool faulty_boolean_generic(unsigned int no, unsigned int input);

extern bool boolean_2(unsigned int no, unsigned int input);

extern float fep[INPUTX];
/*extern bool n[9];
extern bool sa0[9];
extern bool sa1[9];
extern bool on[9];*/
extern unsigned int neg;
extern unsigned int sa0;
extern unsigned int sa1;
extern unsigned int on;

extern unsigned int permutation_no;
extern unsigned int permutation_index;
extern unsigned char permute_order[40321][10];
extern void permute(int level, char *permuted, bool *used, char *original);


#endif // BOOLEAN_H
