#include "fourier.h"

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
