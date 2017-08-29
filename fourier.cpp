#include<math.h>
#include "fourier.h"
#include "boolean.h"

char hdm8[8][8]={
     1,  1,  1,  1,  1,  1,  1,  1,
     1, -1,  1, -1,  1, -1,  1, -1,
     1, 1, -1, -1,  1 , 1, -1, -1,
     1, -1, -1,  1,  1, -1, -1,  1,
     1,  1,  1,  1, -1, -1, -1, -1,
     1, -1,  1, -1, -1,  1, -1,  1,
     1,  1, -1, -1, -1, -1,  1,  1,
     1, -1, -1,  1, -1,  1,  1, -1,


};
char hdm16[16][16], hdm32[32][32], hdm64[64][64], hdm128[128][128], hdm256[256][256];
double fourier_coef[65535];
double influence[15];
double coupling[15];
double total_influence;
double C=6.278;

double max(double a, double b)
{
    if(a>=b) return a;
    else return b;
}

void initialize_fourier()
{
    int i,j;
    /*  Compute Walsh/Hadamard Matrices   *//////////////////////////////////////////////
    for(i=0;i<16;i++){
     for(j=0;j<16;j++){
       hdm16[i][j]=((i>7&&j>7)?(-1):(1))*hdm8[i%8][j%8];
     }
    }
    for(i=0;i<32;i++){
     for(j=0;j<32;j++){
       hdm32[i][j]=((i>15&&j>15)?(-1):(1))*hdm16[i%16][j%16];
     }
    }
    for(i=0;i<64;i++){
     for(j=0;j<64;j++){
       hdm64[i][j]=((i>31&&j>31)?(-1):(1))*hdm32[i%32][j%32];
     }
    }
    for(i=0;i<128;i++){
     for(j=0;j<128;j++){
       hdm128[i][j]=((i>63&&j>63)?(-1):(1))*hdm64[i%64][j%64];
     }
    }
    for(i=0;i<256;i++){
     for(j=0;j<256;j++){
       hdm256[i][j]=((i>127&&j>127)?(-1):(1))*hdm128[i%128][j%128];
     }
    }

}

void calculate_fourier(unsigned int x,bool (*bool_func)(unsigned int,unsigned int)){
    unsigned i,j,k,m;
    double a;
    char fn;
    double nx;

    total_influence=0;

    for(i=0; i<pow(2,INPUT); i++) {
       a=0;
       for(j=0; j<pow(2,INPUT); j++){

           nx=1.0;
           //if(((j^i)|i)==j)
           for(m=0;m<INPUT;m++){
                if(((i>>m)&0x1)==1 && ((j>>m)&0x1)==1) {

                    nx*=(-1.0); //*probability[m];

                } //else nx*=(1-probability[m]);
           }

        fn=bool_func(x,j);
        if(fn==1)fn=-1; else if(fn==0) fn=1;
        //a+=(hdm256[i][j]*fn);
        a+=(nx*fn);
       }
       //if(a!=0) printf(KRED); else printf(KWHT);//
       //if(i==1||i==2||i==4||i==8||i==16) printf(KGRN); else printf(KWHT);//

       //if(a==0)printf("----- "); else printf("%.3f ",(a/16.0));

       fourier_coef[i]=a/pow(2,INPUT);    //pow(2,INPUT);

       for(k=0;k<INPUT;k++)
          if(i==(unsigned)pow(2,k)) {
              influence[k]=fabs(fourier_coef[i]);
              total_influence+=influence[k];
          }

     }

}

double spectral_entropy(){

    int i;
    double h=0;

    for(i=0; i<pow(2,INPUT); i++){

      h+=(pow(fourier_coef[i],2)*log2(1.0/pow(fourier_coef[i],2)));


    }

    return h;

}
double calculate_influence(unsigned int j)
{
 int i;
 double nx=0;
 for(i=0; i<pow(2,INPUT);i++){
     if(((i>>j)&0x1)==1) nx+=pow(fourier_coef[i],2.0);

 }
 return nx;//pow(2,INPUT);
}

double calculate_influence2(unsigned int j)
{ /* using fourier coefficients */
 int i;
 double nx=0;
 for(i=0; i<pow(2,INPUT);i++){
     if(((i>>j)&0x1)==1) nx+=pow(fourier_coef[i],2.0);

 }
 return nx;//pow(2,INPUT);
}


double calculate_coupling(unsigned int i, unsigned int j,bool (*bool_func)(unsigned int,unsigned int))
{
    int k;
    double x=0,x1=0,x2=0;

        for(k=0; k<pow(2,INPUT);k++){

            if(bool_func(0,k)!=bool_func(0,k^(0x1<<j))) x1++;
            if(bool_func(0,k)!=bool_func(0,k^(0x1<<i))) x2++;

            if(bool_func(0,k)!=bool_func(0,k^(0x1<<j)) && bool_func(0,k)!=bool_func(0,k^(0x1<<i)) ) x++;

        }

    return x/max(x1,x2);
    //  return x/pow(2,INPUT);
}

double calculate_coupling2(unsigned int i, unsigned int j)
{ /* using fourier coefficients */
    int k,k2;
    double nx=0;
    for(k=0; k<pow(2,INPUT);k++){

       //for(k2=0; k2<pow(2,INPUT);k2++)
        if(((k>>i)&0x1)==1 && ((k>>j)&0x1)==1) nx+=pow(fourier_coef[k],2.0);

    }
    return nx/max(calculate_influence2(i),calculate_influence2(j));//pow(2,INPUT);
   //return nx/pow(2,INPUT);

}

double total_coupling(unsigned int i){

    double tot=0;

   // return calculate_coupling2(i,0);

    for(int k=0; k<INPUT; k++){
        if(k==i) continue;
        tot+=calculate_coupling2(k,i);
    }
    return tot;
}

/* Degree of f */
unsigned int degree()
{
    int i,j;
    int nx=0, max=0;
    for(i=0; i<pow(2,INPUT);i++){

        for(j=0, nx=0; j<INPUT; j++)
            if(((i>>j)&0x1)==1) nx++;

        if(nx>max) max=nx;

    }
    return max;
}
/* Degree-k weight of f */
double levelWeight(unsigned k)
{
    int i,j;
    int nx=0;
    double tot=0;
    for(i=0; i<pow(2,INPUT);i++){

        for(j=0, nx=0; j<INPUT; j++)
            if(((i>>j)&0x1)==1) nx++;

        if(nx==k) tot+=pow(fourier_coef[i],2);

    }
    return tot;
}
double Stab(double r)
{
    int i,j;
    int nx=0;
    double tot=0;

    /*for(i=0; i<pow(2,INPUT);i++){

        for(j=0, nx=0; j<INPUT; j++)
            if(((i>>j)&0x1)==1) nx++;

        tot+=(pow(r,nx)*pow(fourier_coef[i],2));

    }*/
    for(i=0; i<=INPUT;i++){

            tot+=(pow((1.0-2.0*r),i)*levelWeight(i));

        }
    return tot;
}
double NoiseSensitivity(double r)
{
    return 0.5-0.5*Stab(r);
}

