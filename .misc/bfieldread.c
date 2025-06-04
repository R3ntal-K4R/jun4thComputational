#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG
#define DEBUG1
#define DEBUGBOUNDS


int main(){
    int numlines,offset;
    double x,y,z;
    double delta =.1, deltax = .1, deltay = .1, deltaz = .1; //in meters
	FILE * bfile;
	bfile = fopen("Bfield_n=100.dat","r");
    FILE * wcfile = popen("wc -l Bfield_n=100.dat","r");
    fscanf(wcfile,"%d",&numlines);
    printf("Numlines=%d\n",numlines);
    int dim = (int)cbrt(numlines);
    printf("dim=%d\n",dim);
    printf("Attempting to allocate with size(bytes):%lu\n",numlines*3*sizeof(double));

    
    //Main Bfield array
#ifdef DEBUG
    printf("Allocating array...\n");
    fflush(stdout);
#endif
    double (* BFx)[dim][dim];
    BFx = malloc(dim * sizeof *BFx);
    
    double (* BFy)[dim][dim];
    BFy = malloc(dim * sizeof *BFy);
    
    double (* BFz)[dim][dim];
    BFz = malloc(dim * sizeof *BFz);
    
//    exit(0);
//    double Bfieldx[dim][dim][dim];
//    double Bfieldy[dim][dim][dim];
//    double Bfieldz[dim][dim][dim];
 
    //CHECK BOUNDS
#ifdef DEBUGBOUNDS
    for (int i=0; i<dim; i++){
        printf("BOUNDS CHECK: i=%d",i);
    printf("%.4e %12.4e %12.4e\n",BFx[i][i][i],BFy[i][i][i],BFz[i][i][i]);
    }
#endif
    
#ifdef DEBUG
    printf("Allocated Arrays! Now reading file..\n");
#endif
	//Read all lines
    int i,j,k;
    double Bx,By,Bz;
    fscanf(bfile,"%d",&offset);
    printf("Offset=%d\n",offset);
    rewind(bfile);
    while(fscanf(bfile,"%d %d  %d %lf %lf %lf",&i,&j,&k,&Bx,&By,&Bz) != EOF){
#ifdef DEBUG1
           printf("%d %d  %d %12.4e %12.4e %12.4e\n",i,j,k,Bx,By,Bz);
#endif
        i-=offset;
        j-=offset;
        k-=offset;
        
#ifdef DEBUG1
        printf("x");
        fflush(stdout);
#endif
        BFx[i][j][k] = Bx;
#ifdef DEBUG1
        printf("x!");
        fflush(stdout);
#endif

#ifdef DEBUG1
          printf("y");
        fflush(stdout);
#endif
        BFy[i][j][k] = By;
#ifdef DEBUG1
        printf("y!");
        fflush(stdout);
#endif
#ifdef DEBUG1
          printf("z");
        fflush(stdout);
#endif
        BFz[i][j][k] = Bz;
#ifdef DEBUG1
        printf("z!");
        fflush(stdout);
#endif
        
#ifdef DEBUG1
        printf("IN ARRAY: %d %d  %d %12.4e %12.4e %12.4e\n",i,j,k,BFx[i][j][k],BFy[i][j][k],BFz[i][j][k]);
#endif

        
    }
    printf("Allocation and initialization successful!\n");
    printf("sizeof(BFx[0],BFx[0][0],BFx[0][0][0]), %12.4lu %12.4lu %12.4lu\n",sizeof(BFx[0]),sizeof(BFx[0][0]),sizeof(BFx[0][0][0]));
     printf("sizeof(BFy[0],BFy[0][0],BFy[0][0][0]), %12.4lu %12.4lu %12.4lu\n",sizeof(BFy[0]),sizeof(BFy[0][0]),sizeof(BFy[0][0][0]));
     printf("sizeof(BFz[0],BFz[0][0],BFz[0][0][0]), %12.4lu %12.4lu %12.4lu\n",sizeof(BFz[0]),sizeof(BFz[0][0]),sizeof(BFz[0][0][0]));
    while(1){
    printf("Enter coords x,y,z to poll array\n");
        scanf("%lf %lf %lf",&x,&y,&z);
//        i = (int)fmod(x,delta);
//        j = (int)fmod(y,delta);
//        k = (int)fmod(z,delta);
        i = (int)x/delta;
        j = (int)y/delta;
        k = (int)z/delta;
        
        deltax = fmod(x,delta);
        deltay = fmod(y,delta);
        deltaz = fmod(z,delta);
        
        //Find Bfield with interpolation
//        deltax = x - i*delta;
//        deltay = y - j*delta;
//        deltaz = z - k*delta;
        
        
        
        //Bx(x,y,z)
        Bx = BFx[i][j][k] + (deltax/delta)*(BFx[i+1][j][k]-BFx[i][j][k]) + (deltay/delta)*(BFx[i][j+1][k]-BFx[i][j][k]) + (deltaz/delta)*(BFx[i][j][k+1]-BFx[i][j][k]);
        
        By = BFy[i][j][k] + (deltax/delta)*(BFy[i+1][j][k]-BFy[i][j][k]) + (deltay/delta)*(BFy[i][j+1][k]-BFy[i][j][k]) + (deltaz/delta)*(BFy[i][j][k+1]-BFy[i][j][k]);
        
        Bz = BFz[i][j][k] + (deltax/delta)*(BFz[i+1][j][k]-BFz[i][j][k]) + (deltay/delta)*(BFz[i][j+1][k]-BFz[i][j][k]) + (deltaz/delta)*(BFz[i][j][k+1]-BFz[i][j][k]);
        
#ifdef DEBUG
        printf("    x          y               z        i,j, k,      delta,   deltax       deltay,      deltaz\n");
        printf("%12.4e %12.4e %12.4e   %d %d %d %12.4e %12.4e %12.4e %12.4e",x,y,z,i,j,k, delta, deltax, deltay, deltaz);
#endif
        
        
        printf("\nBx, By, Bz: %12.4e %12.4e %12.4e\n",Bx,By,Bz);
        
        
        
        
        
        
        
        
        //    scanf("%d %d %d",&i,&j,&k);
//        i-=offset;
//        j-=offset;
//        k-=offset;
//    printf("%12.4e %12.4e %12.4e\n",BFx[i][j][k],BFy[i][j][k],BFz[i][j][k]);
    }
    free(BFx);
    free(BFy);
    free(BFz);
	return 0;
}
