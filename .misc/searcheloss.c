#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#define RESET   "\033[0m"
#define BOLDRED  "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW   "\033[1m\033[33m"      /* Bold Green */

unsigned long int InitElossFile(char * filename){
    int lines =0;
     FILE * elossfile=fopen(filename, "r");

        printf("File opened:%s\n", filename);//test
     char chr=getc(elossfile);
        while(chr!=EOF)//finds the number of lines of the file to determine the length of the variable arrays
        {
            if(chr=='\n')
            {
                lines=lines+1;
    //            printf("%i\n", lines);
            }

            chr=getc(elossfile);
        }
        fclose(elossfile);

    printf("Program found %d lines in file %s.\n", lines, filename);//test

    return lines-1;
}

void AllocateElossArrays(char * filename,double * energy,double * electronic,double * nuclear,double * total, int lines){
        FILE *elossfile=fopen(filename, "r");
    char chr;
     
#ifdef DEBUG
        printf("File opened:%s\n", filename);//test
#endif
     while((chr=getc(elossfile))!='\n')//read first line
        {}
      //  fscanf(elossfile,"%[^\n]", c); //skips first line
           printf("Eloss File read initially, now allocating eloss arrays\n");//test
    double energy_holder, electronic_holder, nuclear_holder;
        for(int i=0;i<lines;i++)//assigns values to variables
        {
            //fscanf(elossfile, "%lf %lf %lf %lf %lf", &energy_holder, &electronic_holder, &nuclear_holder, &longitude_stragg_holder, &lateral_strag_holder);
            fscanf(elossfile, "%lf %lf %lf", &energy_holder, &electronic_holder, &nuclear_holder);//
            energy[i]=energy_holder;
            electronic[i]=electronic_holder;
            nuclear[i]=nuclear_holder;
            total[i]=electronic[i]+nuclear[i];
    //        printf("scanning line\n");//test
        }
     fclose(elossfile);

    int i;
    while(1){
        printf("Enter index to poll Eloss arrays, -1 to continue:\n");
        scanf("%d",&i);
        if(i == -1){
            return;}
        printf("%lf %lf %lf %lf",energy[i],electronic[i],nuclear[i],total[i]);
    }

    
}
int main(int argc, char ** argv){
     double energy_holder, electronic_holder,nuclear_holder;
     int lines = InitElossFile(argv[1]);
     double * energy=malloc((lines)*sizeof(double));
     double * electronic=malloc((lines)*sizeof(double));
     double *nuclear=malloc((lines)*sizeof(double));
     double * total=malloc((lines)*sizeof(double));
     AllocateElossArrays(argv[1],energy,electronic,nuclear,total,lines);
    
    printf("Yeet\n");
    
    double currE = 1000000.0, currloss;
    int index;
    int currindex=0;
    //Find initial currindex
    
    while(energy[currindex] < currE){
        currindex++;
    }
    currindex--; //Now points to immediate left of currE
    while(1){
        printf("Now currindex = %d, energy[%d] = %lf, currE = %lf, energy[%d] = %lf",currindex,currindex,energy[currindex], currE,currindex+1,energy[currindex+1]);
        printf("Enter loss. currE = %lf, currindex = %d\n", currE, currindex);
        scanf("%lf",&currloss);
        currE -= currloss;
        
        //Update currindex
        while(currE <= energy[currindex]){
            currindex--;
        }
        
    }
    
    
}

