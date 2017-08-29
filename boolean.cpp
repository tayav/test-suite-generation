#include"boolean.h"
#include<stdio.h>
#include<math.h>

unsigned INPUT= 4; //MAX 8


unsigned int permutation_no=0;
unsigned int permutation_index=0;
unsigned char permute_order[40321][10];


int i,j;

void initialize_boolean()
{
   //for(i=0;i<INPUT;i++) {s[i]=i; used[i]=false;}


}

void permute(int level, char *permuted, bool *used, char *original) {
    char p[130]="";
    int i;
    unsigned k;
    if(level == INPUT) { // permutation complete, display
        /*if(permutation_index++==permutation_no){
            for(i=0;i<INPUT;i++)permute_order[i]=permuted[i];
            return;
        }*/

        k=permutation_index;permutation_index++;
        for(i=0;i<INPUT;i++) permute_order[k][i]=permuted[i];
        //for(i=0;i<INPUT;i++) printf("%u",permute_order[k][i]);
        //printf("\n\n");


    } else {
        for(i=0; i<INPUT; i++) { // try to add an unused character
            if(!used[i]) {
                used[i] = true;
                sprintf(p,"%c%s",original[i],permuted);
                permute(level+1,  p, used, original); // find the permutations starting with this string
                used[i] = false;
            }
        }
}
}

unsigned fact(unsigned int n)
{
    unsigned int i, fact=1;
    for (i = 1; i <= n; i++)
        fact = fact * i;

    return fact;
}

bool opr2(bool x1,bool x2, int opr)
{
    if(opr==AND) return x1*x2;
    else if(opr==OR) return x1+x2;
}

void print_binary(unsigned int k)
{   int i;
    for(i=0;i<INPUT;i++) {
            printf("%u",k>>i&0x1);
    }
}

void fprint_binary(FILE *fp, unsigned int k)
{   int i;
    for(i=0;i<INPUT;i++) {
            fprintf(fp,"%u ",k>>i&0x1);
    }
}

float correlation(float x[], float y[])
{
    float ab,a,b,a2,b2,t;
    for(i=0,t=0;i<INPUT;i++)t+=(x[i]*y[i]); ab=t;
    for(i=0,t=0;i<INPUT;i++)t+=(x[i]); a=t;
    for(i=0,t=0;i<INPUT;i++)t+=(y[i]); b=t;
    for(i=0,t=0;i<INPUT;i++)t+=(x[i]*x[i]); a2=t;
    for(i=0,t=0;i<INPUT;i++)t+=(y[i]*y[i]); b2=t;

    return ((float)INPUT*ab-a*b)/sqrt(((float)INPUT*a2-a*a)*((float)INPUT*b2-b*b));

}

bool boolean_2(unsigned int no, unsigned int input)
{
    bool x0,x1,x2,x3;
    x0=(input>>0)&0x1;
    x1=(input>>1)&0x1;
    x2=(input>>2)&0x1;
    x3=(input>>3)&0x1;


    return (x1*x2)*x0*x3;

}

bool boolean_generic(unsigned int no, unsigned int input)
{
    unsigned i,j;
    //unsigned int total_permutations;
    bool b;

    // change order of input bits : factorial of numofinputs
    // AND / OR for numofinputs-1  operators between the input bits
    // negation of itself for numofinputs bits

    // no:  log(2,n!) bit input permutations | (n-1) bit operations | n bit negations

    //total_permutations=ceil(log2(fact(INPUT)));

    //printf("fact=%u   %u\n bits",fact(numofinputs),total_permutations);

    permutation_no=no>>INPUT>>(INPUT-1);
    //printf("permutation no=%u\n",permutation_no);

    //permutation_index=0;
    //permute(0, p, used, s);

    //for(i=0;i<INPUT;i++)printf("%u",permute_order[i]); printf("\n");

    b=((input>>(permute_order[permutation_no][0]-1))&0x1)^((no>>0)&0x1);

    for(i=1;i<INPUT;i++){

        b=opr2(((input>>(permute_order[permutation_no][i]-1))&0x1)^((no>>i)&0x1),b,(no>>(INPUT-1)>>i)&0x1);


        //printf("%u",((input>>i) & 0x01));

    }


    return b;
}

// FEP: fault exposing potential; defined for each input = mcdc test pair
float fep[INPUTX];

unsigned int neg=0x0000;
unsigned int sa0=0xFFFF;
unsigned int sa1=0x0000;
unsigned int on=0x0000;

/*
bool n[9]=  {0,0,0,0,0,0,0,0,0};
bool sa0[9]={1,1,1,1,1,1,1,1,1};
bool sa1[9]={0,0,0,0,0,0,0,0,0};
bool on[9]={0,0,0,0,0,0,0,0,0};
*/

bool faulty_boolean_generic(unsigned int no, unsigned int input)
{
    unsigned i,j;
    unsigned int total_permutations;
    bool b;

    // change order of input bits : factorial of numofinputs
    // AND / OR for numofinputs-1  operators between the input bits
    // negation of itself for numofinputs bits


    // no:  log(2,n!) bit input permutations | (n-1) bit operations | n bit negations

    total_permutations=ceil(log2(fact(INPUT)));

    //printf("fact=%u   %u\n bits",fact(numofinputs),total_permutations);

    permutation_no=no>>INPUT>>(INPUT-1);
    //printf("permutation no=%u\n",permutation_no);

    //permutation_index=0;
    //permute(0, p, used, s);

    b=(((((input>>(permute_order[permutation_no][0]-1))&0x1)^((no>>0)&0x1))^((neg>>0)&0x1))&((sa0>>0)&0x1))|((sa1>>0)&0x1);



    for(i=1;i<INPUT;i++){

        //b=opr2(((((((input>>(permute_order[permutation_no][i]-1))&0x1)^((no>>i)&0x1))^((neg>>i)&0x1)&((sa0>>i)&0x1))|((sa1>>i)&0x1))),b,((no>>(INPUT-1)>>i)&0x1)^((on>>(i-1))&0x1));
        b=opr2(((((((input>>(permute_order[permutation_no][i]-1))&0x1)^((no>>i)&0x1))^((neg>>i)&0x1)))),b,((no>>(INPUT-1)>>i)&0x1)^((on>>(i-1))&0x1));

        //printf("%u",((input>>i) & 0x01));

    }


    return b;
}
