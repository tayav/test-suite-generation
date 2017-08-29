/*  MUTATION BASED OPTIMUM TEST SUITE GENERATION */

#define VNF_ENF_SINGLE_FAULT
#define ORF_SINGLE_FAULT
#define CCF_SINGLE_FAULT
#define CDF_SINGLE_FAULT


//#define VNF_ENF_DOUBLE_FAULT
//#define ORF_DOUBLE_FAULT
//#define SA0_SINGLE_FAULT
//#define SA1_SINGLE_FAULT
//#define SA0_DOUBLE_FAULT
//#define SA1_DOUBLE_FAULT
//#define CCF_DOUBLE_FAULT
//#define CDF_DOUBLE_FAULT
//#define VRF_SINGLE_FAULT
//#define VRF_DOUBLE_FAULT

//#define FIND_MASKING_MCDC_TESTS
//#define FIND_UC_MCDC_TESTS

#include <QCoreApplication>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include"fourier.h"
#include"boolean.h"
#define FILENAME1 "amatrix.csv"
#define FILENAME2 "coverage.mod"
#define FILENAME3 "coverage.final.sol"
#define FILENAME4 "coverage2.final.sol"
#define FILENAME5 "influences.txt"
#define FILENAME6 "faultdetection.txt"
#define FILENAME7 "faultdetection2.txt"
unsigned TCASNO=0;

bool fx[65535];
bool mcdc_uc[65535];
bool mcdc_masking[65535];
bool x[16];
bool y[100];
unsigned mcdc_test[30][2];

char f[21][100][9];
int oprnum[21]={3,22,35,35,4,19,27,20,16,9,14,15,16,12,11,17,24,10,10,8,7};
int inpnum[21]={4,8 ,9, 13,5,12,11,11,8,7,13,13,14,12,8, 9  ,13,11,11,8,7};
unsigned char amatrix[500][16384];
char s[30],s2[30];
int totalkilled=0,xi,yi;
unsigned tests[50], ordered[50], prio[50];
unsigned totaltests;

extern void loadfdata();

bool opr(uchar i1type, uchar i1, bool i1n, uchar i2type, uchar i2, bool i2n, uchar operation, uchar o, bool on){

    bool a, b, z;

    if(i1type == 'x') a=x[i1]^i1n; else a=y[i1]^i1n;
    if(i2type == 'x') b=x[i2]^i2n; else b=y[i2]^i2n;

    switch(operation){
        case '*': z=a*b; break;
        case '+': z=a+b; break;
        default: printf("operator error: "); exit(1); break;
    }

    y[o]=(z^on);

    return y[o];
}


unsigned negfault[512];
unsigned sa1fault[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
unsigned sa0fault[16]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
unsigned ccf1=0, ccf2=0, cdf1=0, cdf2=0;
unsigned ccf11=0, ccf22=0, cdf11=0, cdf22=0;

bool boolean_tcas(unsigned int n, unsigned int input)
{
    int i;
    bool z;
    unsigned findex=0;

    for(i=0;i<16;i++) x[i]=((((input>>i)&0x1) & sa0fault[i]) | sa1fault[i]);

    x[ccf1]=x[ccf1] & x[ccf2];
    x[ccf11]=x[ccf11] & x[ccf22];

    x[cdf1]=x[cdf1] | x[cdf2];
    x[cdf11]=x[cdf11] | x[cdf22];


    for(i=0;i<oprnum[n];i++){
         //printf("%c %d %d %c %d %d %c %d %d \n",f[n][i][0],f[n][i][1]-48,f[n][i][2]-48,f[n][i][3],f[n][i][4]-48,f[n][i][5]-48,f[n][i][6],f[n][i][7]-48,f[n][i][8]-48);
         z=opr(f[n][i][0],f[n][i][1]-0,(f[n][i][2]-0)^negfault[findex],f[n][i][3],f[n][i][4]-0,(f[n][i][5]-0)^negfault[findex+1],f[n][i][6],f[n][i][7]-0,(f[n][i][8]-0)^negfault[findex+2]);
         findex+=3;
    }

    return z;
}

bool masking_mcdc_check()
{
    bool chk;
    unsigned jj;

    for(jj=0, chk=false; jj<65535; jj++)if(mcdc_masking[jj]){


            if( (boolean_tcas(TCASNO,jj)!= fx[jj]) )
             {
                chk=true;
                return true;
             }

    }

    return chk;
}

void order(unsigned tests[], unsigned prio[], unsigned total)
{
 unsigned int i,j,temp;

 for(i=0; i<total-1; i++)
     for(j=i+1; j<total; j++){

        if(prio[i]<prio[j]){

            temp=prio[j]; prio[j]=prio[i]; prio[i]=temp;
            temp=tests[j]; tests[j]=tests[i]; tests[i]=temp;

        }

     }


}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    unsigned int i,j,k;
    unsigned int aa, bb, cc, dd;
    unsigned int xx,yy;
    unsigned int ii,jj1, jj2, jj;
    unsigned iii,iiii;
    bool tmp, tmp2, chk, chk2;
    char tmpx, tmpx2;
    unsigned mcdc_detected=0, mcdc_notdetected=0;
    int mutantno=1;

    time_t start,end;

    FILE *fp=NULL, *mp=NULL, *ip=NULL, *tp=NULL, *tp2=NULL;

    if(argc<2) {printf("Missing parameter: TCAS# (1-20) must be specified.\n");exit(2);}
    else TCASNO=atoi(argv[1]);

    fp=fopen(FILENAME1,"w");
    fprintf(fp,"ROWX,COLX,VALX\n");

    time(&start);


 //int tcasx=1;
 //for(tcasx=4;tcasx<=4   ;tcasx++){

     //TCASNO=tcasx;

    INPUT=inpnum[TCASNO];

    for(i=0;i<512;i++) negfault[i]=false;
    for(i=0;i<65535;i++) mcdc_masking[i]=mcdc_uc[i]=false;

    loadfdata();

    printf("Utility for Mutation-based Optimum Structural Coverage\n------------------------------------------------------\n");
    printf("TCAS-II Specifications\n");

    initialize_fourier();

    printf("TCAS # %d\nTruth vector:\n",TCASNO);//for(i=0;i<pow(2,inpnum[TCASNO]);i++) printf("%u ",boolean_tcas(TCASNO,i));

    printf("\n# of inputs (n):%u\n",inpnum[TCASNO]);
    printf("# of operations:%u\n",oprnum[TCASNO]);



    for(i=0;i<pow(2,inpnum[TCASNO]);i++) fx[i]=boolean_tcas(TCASNO,i);

    calculate_fourier(TCASNO,boolean_tcas);

    for(i=0;i<inpnum[TCASNO];i++) influence[i]=calculate_influence(i);

    for(i=0;i<inpnum[TCASNO];i++) coupling[i]=total_coupling(i);
    //fprintf(stdout,"\nFour:\t");
    //for(i=0;i<pow(2,inpnum[TCASNO]);i++) fprintf(stdout,"%f ",fourier_coef[i]);
    fprintf(stdout,"\nInf:\t");
    for(i=0;i<inpnum[TCASNO];i++) fprintf(stdout,"%f ",influence[i]); fprintf(stdout,"   I[f]:%.4f ",total_influence);

    ip=fopen(FILENAME5,"a+");
    fprintf(ip,"%u\t",TCASNO);
    for(i=0;i<inpnum[TCASNO];i++) fprintf(ip,"%f\t",influence[i]); fprintf(ip,"\t%.4f\n",total_influence);
    fclose(ip);

    fprintf(stdout,"\nCoup:\t");
    //for(i=0;i<inpnum[TCASNO];i++) fprintf(stdout,"%f\t",coupling[i]);
    //fprintf(stdout,"\nStability: %f\n",Stab(0.1));

    /* Find out Masking MCDC test pairs *//////////////////////////////////////////////////////
    #ifdef FIND_MASKING_MCDC_TESTS
    printf("Finding Masking MCDC test cases "); fflush(stdout);
    for(i=0;i<inpnum[TCASNO];i++){ printf(".");fflush(stdout);
    for(xx=0;xx<pow(2,inpnum[TCASNO]);xx++)for(yy=0;yy<pow(2,inpnum[TCASNO]);yy++){
        //if(yy>=xx) continue;
        if((((xx>>i)&0x1)^((yy>>i)&0x1)) && (boolean_tcas(TCASNO,xx)^boolean_tcas(TCASNO,yy)) &&
                (boolean_tcas(TCASNO,xx^(1<<i))^boolean_tcas(TCASNO,xx))  &&
                (boolean_tcas(TCASNO,yy^(1<<i))^boolean_tcas(TCASNO,yy))
                  && xx==(yy^(1<<i))
                // bits do not change in x and y? requirement of masking MCDC? What is the effect of this?
                ){

            mcdc_test[i][0]=xx;
            mcdc_test[i][1]=yy;

            mcdc_masking[xx]=true;
            mcdc_masking[yy]=true;

            xx=yy=65535;

            //printf("\n");
            //printf("i=%u\t",i); print_binary(x);printf("\t");print_binary(y);
            //printf("\t(%u\t%u)\n",x,y);
            break;
        }
    }    }
    printf("OK\n"); fflush(stdout);
    printf("MCDC pairs:\t"); for(jj=0;jj<inpnum[TCASNO];jj++) printf("(%u,%u)\t",mcdc_test[jj][0],mcdc_test[jj][1]);
    printf("\nMasking MCDC tests:\t"); for(jj=0,i=0;jj<pow(2,inpnum[TCASNO]);jj++) if(mcdc_masking[jj]) { printf("%u ",jj); i++;}
    printf("\n");
    for(jj=0;jj<pow(2,inpnum[TCASNO]);jj++) if(mcdc_masking[jj]) { fprint_binary(stdout,jj); printf(" :%u  (%u)\n",fx[jj],jj);}
    printf("\n# of Masking MCDC tests: %u\t",i);
    #endif

    #ifdef FIND_UC_MCDC_TESTS
    printf("\nFinding Unique Cause MCDC test cases "); fflush(stdout);
    for(i=0;i<inpnum[TCASNO];i++){ printf(".");fflush(stdout);
    for(xx=0;xx<pow(2,inpnum[TCASNO]);xx++)for(yy=0;yy<pow(2,inpnum[TCASNO]);yy++){
        //if(yy>=xx) continue;
        if((((xx>>i)&0x1)^((yy>>i)&0x1)) && (boolean_tcas(TCASNO,xx)^boolean_tcas(TCASNO,yy)) &&
                (boolean_tcas(TCASNO,xx^(1<<i))^boolean_tcas(TCASNO,xx))  &&
                (boolean_tcas(TCASNO,yy^(1<<i))^boolean_tcas(TCASNO,yy))
                  && xx==(yy^(1<<i))
                // bits do not change in x and y? requirement of Unique-cause MCDC? What is the effect of this?
                ){

            //if(xx==0 && yy==0) continue;

            //mcdc_test[i][0]=xx;
            //mcdc_test[i][1]=yy;

            mcdc_uc[xx]=true;
            mcdc_uc[yy]=true;

            //xx=yy=65535;

            //printf("\n");
            //printf("i=%u\t",i); print_binary(x);printf("\t");print_binary(y);
            //printf("\t(%u\t%u)\n",x,y);
            break;
        }
    }    }
    printf("OK\n"); fflush(stdout);
    //printf("UC MCDC pairs:\t"); for(jj=0;jj<inpnum[TCASNO];jj++) printf("(%u,%u)\t",mcdc_test[jj][0],mcdc_test[jj][1]);
    printf("\nUC MCDC tests:\t"); for(jj=0,i=0;jj<pow(2,inpnum[TCASNO]);jj++) if(mcdc_uc[jj]) { printf("%u ",jj); i++;}    
    printf("\n# of UC MCDC tests: %u\t",i);

    #endif
    /* END of Find out MCDC test pairs *///////////////////////////////////////////////

    /* Generate Mutants *//////////////////////////////////////////////////////////////

    int eliminatedmutant=0;
    printf("\nGenerating mutants "); fflush(stdout);
#ifdef VNF_ENF_SINGLE_FAULT
    printf(".");fflush(stdout);
    for(i=0;i<oprnum[TCASNO]*3;i++){

        negfault[i]=1;
        //printf("\n");
        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u & ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
        negfault[i]=0; if(chk)mutantno++; else eliminatedmutant++;

    }
#endif

#ifdef VNF_ENF_DOUBLE_FAULT
    printf(".");fflush(stdout);
    for(i=0;i<oprnum[TCASNO]*3;i++)for(ii=0;ii<oprnum[TCASNO]*3;ii++){
        if(i==ii) continue;

        negfault[i]=negfault[ii]=1;
        //printf("\n");
        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
        negfault[i]=negfault[ii]=0; if(chk)mutantno++; else eliminatedmutant++;

    }
#endif



#ifdef VRF_SINGLE_FAULT
    printf(".");fflush(stdout);
    unsigned ri1,ri2,ci1,ci2;
    for(aa=0;aa<oprnum[TCASNO]*2;aa++)for(bb=0;bb<oprnum[TCASNO]*2;bb++){

        if(aa==bb) continue;

        if(aa>=oprnum[TCASNO]) { ri1=aa-oprnum[TCASNO]; ci1=3;} else { ri1=aa; ci1=0;}
        if(bb>=oprnum[TCASNO]) { ri2=bb-oprnum[TCASNO]; ci2=3;} else { ri2=bb; ci2=0;}

        if(f[TCASNO][ri1][ci1]!='x' || f[TCASNO][ri2][ci2]!='x') continue;

        tmpx=f[TCASNO][ri1][ci1+1];
        f[TCASNO][ri1][ci1+1]=f[TCASNO][ri2][ci2+1];
        f[TCASNO][ri2][ci2+1]=tmpx;

        //printf("%u %u - %u %u - %u %u - %u %u %u\n",ri1,ci1,ri2,ci2,f[TCASNO][ri1][ci1+1],f[TCASNO][ri2][ci2+1],aa,bb,tmpx);

        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
        if(chk)mutantno++;else eliminatedmutant++;

        tmpx=f[TCASNO][ri1][ci1+1];
        f[TCASNO][ri1][ci1+1]=f[TCASNO][ri2][ci2+1];
        f[TCASNO][ri2][ci2+1]=tmpx;

    }
#endif
#ifdef VRF_DOUBLE_FAULT
    printf(".");fflush(stdout);
    unsigned ri1,ri2,ci1,ci2;
    for(aa=0;aa<oprnum[TCASNO]*2;aa++)for(bb=0;bb<oprnum[TCASNO]*2;bb++)for(cc=0;cc<oprnum[TCASNO]*2;cc++)for(dd=0;dd<oprnum[TCASNO]*2;dd++){

        if(aa==bb) continue;
        if(cc==dd) continue;

        ri1=(int)floor((float)aa/2);
        ri2=(int)floor((float)bb/2);
        if((aa%2)==0) ci1=0; else ci1=3;
        if((bb%2)==0) ci2=0; else ci2=3;
        if(f[TCASNO][ri1][ci1]!='x' || f[TCASNO][ri2][ci2]!='x') continue;
        tmp=f[TCASNO][ri1][ci1+1];
        f[TCASNO][ri1][ci1+1]=f[TCASNO][ri2][ci2+1];
        f[TCASNO][ri2][ci2+1]=tmp;


        ri1=(int)floor((float)cc/2);
        ri2=(int)floor((float)dd/2);
        if((cc%2)==0) ci1=0; else ci1=3;
        if((dd%2)==0) ci2=0; else ci2=3;
        if(f[TCASNO][ri1][ci1]!='x' || f[TCASNO][ri2][ci2]!='x') continue;
        tmp=f[TCASNO][ri1][ci1+1];
        f[TCASNO][ri1][ci1+1]=f[TCASNO][ri2][ci2+1];
        f[TCASNO][ri2][ci2+1]=tmp;


        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
        if(chk)mutantno++;else eliminatedmutant++;

        ri1=(int)floor((float)aa/2);
        ri2=(int)floor((float)bb/2);
        if((aa%2)==0) ci1=0; else ci1=3;
        if((bb%2)==0) ci2=0; else ci2=3;
        if(f[TCASNO][ri1][ci1]!='x' || f[TCASNO][ri2][ci2]!='x') continue;
        tmp=f[TCASNO][ri1][ci1+1];
        f[TCASNO][ri1][ci1+1]=f[TCASNO][ri2][ci2+1];
        f[TCASNO][ri2][ci2+1]=tmp;

        ri1=(int)floor((float)cc/2);
        ri2=(int)floor((float)dd/2);
        if((cc%2)==0) ci1=0; else ci1=3;
        if((dd%2)==0) ci2=0; else ci2=3;
        if(f[TCASNO][ri1][ci1]!='x' || f[TCASNO][ri2][ci2]!='x') continue;
        tmp=f[TCASNO][ri1][ci1+1];
        f[TCASNO][ri1][ci1+1]=f[TCASNO][ri2][ci2+1];
        f[TCASNO][ri2][ci2+1]=tmp;

    }
#endif



#ifdef SA0_SINGLE_FAULT
    printf(".");fflush(stdout);
    for(i=0;i<inpnum[TCASNO];i++){

        sa0fault[i]=0;        
        //printf("\n");
        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
        sa0fault[i]=1; if(chk)mutantno++; else eliminatedmutant++;
    }
#endif

#ifdef SA0_DOUBLE_FAULT
    printf(".");fflush(stdout);
    for(i=0;i<inpnum[TCASNO];i++)for(ii=0;ii<inpnum[TCASNO];ii++){
        if(i==ii) continue;
        sa0fault[i]=sa0fault[ii]=0;
        //printf("\n");
        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
        sa0fault[i]=sa0fault[ii]=1; if(chk)mutantno++; else eliminatedmutant++;
    }
#endif
#ifdef SA1_SINGLE_FAULT
    printf(".");fflush(stdout);
    for(i=0;i<inpnum[TCASNO];i++){

        sa1fault[i]=1;

        //printf("\n");
        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
        sa1fault[i]=0; if(chk)mutantno++; else eliminatedmutant++;
    }
#endif
#ifdef SA1_DOUBLE_FAULT
    printf(".");fflush(stdout);
    for(i=0;i<inpnum[TCASNO];i++)for(ii=0;ii<inpnum[TCASNO];ii++){
        if(i==ii) continue;
        sa1fault[i]=sa1fault[ii]=1;
        //printf("\n");
        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
        sa1fault[i]=sa1fault[ii]=0; if(chk)mutantno++; else eliminatedmutant++;
    }
#endif

#ifdef CCF_SINGLE_FAULT
    printf(".");fflush(stdout);
    for(i=0;i<inpnum[TCASNO];i++)for(ii=0;ii<inpnum[TCASNO];ii++){
        if(i==ii) continue;
        ccf1=i; ccf2=ii;
        //printf("\n");
        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
         if(chk)mutantno++; else eliminatedmutant++;
    }
    ccf1=ccf2=0;
#endif

#ifdef CCF_DOUBLE_FAULT
    printf(".");fflush(stdout);
    for(i=0;i<inpnum[TCASNO];i++)for(ii=0;ii<inpnum[TCASNO];ii++)for(iii=0;iii<inpnum[TCASNO];iii++)for(iiii=0;iiii<inpnum[TCASNO];iiii++){
        if(i==ii || iii==iiii) continue;
        ccf1=i; ccf2=ii; ccf11=iii; ccf22=iiii;
        //printf("\n");
        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
         if(chk)mutantno++; else eliminatedmutant++;
    }
    ccf1=ccf2=ccf11=ccf22=0;
#endif


#ifdef CDF_SINGLE_FAULT
    printf(".");fflush(stdout);
    for(i=0;i<inpnum[TCASNO];i++)for(ii=0;ii<inpnum[TCASNO];ii++){
        if(i==ii) continue;
        cdf1=i; cdf2=ii;
        //printf("\n");
        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
         if(chk)mutantno++; else eliminatedmutant++;
    }
    cdf1=cdf2=0;
#endif

#ifdef CDF_DOUBLE_FAULT
    printf(".");fflush(stdout);
    for(i=0;i<inpnum[TCASNO];i++)for(ii=0;ii<inpnum[TCASNO];ii++)for(iii=0;iii<inpnum[TCASNO];iii++)for(iiii=0;iiii<inpnum[TCASNO];iiii++){
        if(i==ii || iii==iiii) continue;
        cdf1=i; cdf2=ii;
        cdf11=iii; cdf22=iiii;
        //printf("\n");
        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
         if(chk)mutantno++; else eliminatedmutant++;
    }
    cdf1=cdf2=cdf11=cdf22=0;
#endif

#ifdef ORF_SINGLE_FAULT
    printf(".");fflush(stdout);
    for(i=0;i<oprnum[TCASNO];i++){

        tmpx=f[TCASNO][i][6];
        if(f[TCASNO][i][6]=='*') f[TCASNO][i][6]='+'; else f[TCASNO][i][6]='*';
        //printf("\n");
        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
        f[TCASNO][i][6]=tmpx; if(chk)mutantno++; else eliminatedmutant++;
    }
#endif

#ifdef ORF_DOUBLE_FAULT
    printf(".");fflush(stdout);
    for(i=0;i<oprnum[TCASNO];i++)for(ii=0;ii<oprnum[TCASNO];ii++){
        if(i==ii) continue;

        tmpx=f[TCASNO][i][6];tmpx2=f[TCASNO][ii][6];
        if(f[TCASNO][i][6]=='*') f[TCASNO][i][6]='+'; else f[TCASNO][i][6]='*';
        if(f[TCASNO][ii][6]=='*') f[TCASNO][ii][6]='+'; else f[TCASNO][ii][6]='*';
        //printf("\n");
        for(chk=false, j=0;j<pow(2,inpnum[TCASNO]);j++) {
            if(boolean_tcas(TCASNO,j)^fx[j]) { chk=true;  break; }
        };
        if(chk)for(j=0;j<pow(2,inpnum[TCASNO]);j++) {
            //printf("%u ",boolean_tcas(TCASNO,j)^fx[j]);
            fprintf(fp,"%u,%u,%u\n",mutantno,j+1,boolean_tcas(TCASNO,j)^fx[j]);
        };
        if(chk){ if(masking_mcdc_check()) mcdc_detected++; else mcdc_notdetected++;}
        f[TCASNO][i][6]=tmpx;f[TCASNO][ii][6]=tmpx2; if(chk)mutantno++; else eliminatedmutant++;
    }
#endif

    printf("OK\t %u mutants. %u eliminated. \n",mutantno-1, eliminatedmutant); fflush(stdout);
    /* End of Generate Mutants  *//////////////////////////////////////////////////////
    printf("\n");
    printf("MCDC:   Detected=%u,   Not detected=%u\n",mcdc_detected,mcdc_notdetected);
    printf("MCDC Coverage=%.3f\n",(float)mcdc_detected/((float)mutantno-1));

//}

    mp=fopen(FILENAME2,"w");
    fprintf(mp,"param inputnum, integer, >0;\nparam mutantnum, integer, >0;\n");
    fprintf(mp,"set INPUT := 1..inputnum;\nset MUTANT  := 1..mutantnum ;\nset A dimen 2;\n");
    fprintf(mp,"param a {i in MUTANT, j in INPUT};\ntable data IN \"CSV\" \"amatrix.csv\": A<-[ROWX,COLX],a~VALX;\nvar x {i in INPUT} >=0, binary;\nvar y {i in MUTANT} >=0, <=inputnum-1, integer;\n");

    //fprintf(mp,"minimize z: 0.9*sum {i in INPUT} (  (1/sum {r in MUTANT} (a[r,i])) * x[i]) + 0.1*sum {r in MUTANT}( sum {i in INPUT}(a[r,i]*x[i])) ;\n");
    //fprintf(mp,"s.t. const{r in MUTANT}: sum{i in INPUT: a[r,i]==1} (x[i]) >= 1;\n");

    fprintf(mp,"minimize z: 0.9999*sum {i in INPUT} (  (1/(sum {r in MUTANT} (a[r,i])+1)) * x[i]) + (1-0.9999)*sum {r in MUTANT}(y[r]) ;\n");
    fprintf(mp,"s.t. const{r in MUTANT}: sum{i in INPUT: a[r,i]==1} (x[i]) >= inputnum-y[r];\n");

    //fprintf(mp,"minimize z: sum {i in INPUT} ( x[i]) ;\n");
    //fprintf(mp,"s.t. const{r in MUTANT}: sum{i in INPUT: a[r,i]==1} (x[i]) >= 1;\n");

    fprintf(mp,"solve;\n");

    //fprintf(mp,"display {i in INPUT} sum {r in MUTANT} (a[r,i]) ;\n");
    fprintf(mp,"display sum{r in MUTANT} (y[r]) ;\n");
    fprintf(mp,"printf \"# of Test Cases= %%"); fprintf(mp,"d\\n\",  sum{i in 1..inputnum} x[i];\n");
    fprintf(mp,"\n");
    fprintf(mp,"table result{i in INPUT} OUT \"CSV\" \"coverage.final.sol\": i~ROW,x[i]~VALX,sum {r in MUTANT} (a[r,i])~PRIOX;\n");
    fprintf(mp,"table result{i in INPUT} OUT \"CSV\" \"coverage2.final.sol\": i~ROW,sum {r in MUTANT}(a[r,i])~VALX;\n");
    fprintf(mp,"data;\n");
    fprintf(mp,"param inputnum := %u;\n",(unsigned int)pow(2,inpnum[TCASNO]));
    fprintf(mp,"param mutantnum := %u;\n",mutantno-1);
    fprintf(mp,"end;\n");
    fcloseall();

    printf("\nPress any key to start the optimization.");
    getchar();
    time(&end);
    system("glpsol -m coverage.mod");


    printf("\n# of mutants : %u\n",mutantno-1);
    printf("Time=%.3f minutes\nTest Cases:\n",(double)(end-start)/60.0);



    ii=0;

    fp=fopen(FILENAME3,"r");
    while(fgets(s,100,fp)){
        sscanf(s,"%u,%u,%u\n",&i,&j,&k);
        if(j==1) { fprint_binary(stdout,i-1); printf(":%u (%u) --> %u ",boolean_tcas(TCASNO,i-1),i-1,k); printf("\n");}


    };
    fclose(fp);
    fp=fopen(FILENAME3,"r");
    mp=fopen(FILENAME4,"r");
    while(fgets(s,100,fp)){
        sscanf(s,"%u,%u\n",&i,&j);
        if(fgets(s2,100,mp)) sscanf(s2,"%u,%u\n",&xi,&yi);



        if(j==1) {  /*printf("%u(%u), ",i-1,yi); */ totalkilled+=yi;

                    tests[ii]=i-1; prio[ii]=yi; ii++;
        }


    };

    totaltests=ii;

    for(j=0; j<ii; j++) printf("%u(%u),",tests[j],prio[j]);

    //order(tests,prio,ii);
    //printf("\n");
    //for(j=0; j<ii; j++) printf("%u(%u),",tests[j],prio[j]);

    printf("\nTotal Killed by Tests = %d\n",totalkilled);

    /* Calculate APFD */

    fp=fopen(FILENAME1,"r");
    fgets(s,100,fp);
    while(fgets(s,100,fp)){
        sscanf(s,"%u,%u,%u\n",&i,&j,&k);
        amatrix[i-1][j-1]=k;
    }

    tp=fopen(FILENAME6,"w");
    tp2=fopen(FILENAME7,"w");
    unsigned totalTF=0;
    float apfd;

    for(i=0; i<mutantno-1; i++) {
        for(j=0; j<totaltests; j++){

             if(amatrix[i][tests[j]]==1) { totalTF+=(j+1); fprintf(tp,"%u\n",j+1); break; }

        }
    }
    apfd=1.0-((float)totalTF/(float)(totaltests*(mutantno-1)))+1.0/(2.0*totaltests);

    printf("APFD=%.3f\n",apfd);

    order(tests,prio,totaltests);
    printf("\n");
    for(j=0; j<totaltests; j++) printf("%u(%u),",tests[j],prio[j]);
    printf("\n");



    totalTF=0;
    for(i=0; i<mutantno-1; i++) {
        for(j=0; j<totaltests; j++){

             if(amatrix[i][tests[j]]==1) { totalTF+=(j+1); fprintf(tp2,"%u\n",j+1);  break; }

        }
    }
    apfd=1.0-((float)totalTF/(float)(totaltests*(mutantno-1)))+1.0/(2.0*totaltests);

    printf("APFD*=%.3f\n",apfd);

    exit(1);
    return a.exec();
}
