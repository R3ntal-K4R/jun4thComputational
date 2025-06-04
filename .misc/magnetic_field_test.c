#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "aux.h"
#include "mag.h"
#include "eloss.h"
#include "motion.h"

void magnetic_field_test(const char *filename){
    FILE *bfile;
    bfile = fopen(filename, "r");
    char bfieldcmd[100] = "wc -l ";
    strcat(bfieldcmd, filename);
    printf("Using file %s\n", filename);

    int numlines;
    FILE * wcfile = popen(bfieldcmd, "r");
    fscanf(wcfile, "%d", &numlines);

    printf("%s\n", bfieldcmd);
    printf("Bfield Numlines = %d\n", numlines);

    int dim;
    dim = (int)cbrt(numlines);
    printf("Dimensions = %d\n", dim);

    FILE *ptfilein; // input file
    ptfilein = fopen("mag_field_test_input.dat", "r");
    if (ptfilein == NULL) {
        printf("Error: No input file! Make sure a file named 'mag_field_test_input.dat' is present in the current directory\n");
        exit(0);
    }

    double x, y, z;
    double x_plane, y_plane, z_plane;

    fscanf(ptfilein, "%lf %lf %lf", &x, &y, &z);
    fscanf(ptfilein, "%lf %lf %lf", &x_plane, &y_plane, &z_plane);
    fclose(ptfilein);

    double static (* BFx)[dim][dim];
    BFx = malloc(dim * sizeof *BFx);

    double static (* BFy)[dim][dim];
    BFy = malloc(dim * sizeof *BFy);
        
    double static (* BFz)[dim][dim];
    BFz = malloc(dim * sizeof *BFz);

    int i,j,k, offset;
    double currBx,currBy,currBz;
    fscanf(bfile,"%d",&offset);

    rewind(bfile);
     while(fscanf(bfile,"%d %d  %d %lf %lf %lf",&i,&j,&k,&currBx,&currBy,&currBz) != EOF){
        i-=offset;
        j-=offset;
        k-=offset;

        BFx[i][j][k] = currBx;
        BFy[i][j][k] = currBy;
        BFz[i][j][k] = currBz;
     }

    printf("Allocated Arrays! Now creating data files.\n");

    FILE *bfield_test_file = fopen("bfield_test_file.dat", "w");

    if(x_plane == 1){
        printf("Creating data file for a heat map for a y-z plane at x = %lf!\n", x);
        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                double curr_mag_of_bfield = rmag(Bx(BFx,x,i,j),By(BFy,x,i,j),Bz(BFz,x,i,j));
                fprintf(bfield_test_file, "%d %d %lf\n", i, j, curr_mag_of_bfield);
            }
        }
    }
    else if(y_plane == 1){
        printf("Creating data file for a heat map for an x-z plane at y = %lf!\n", y);
        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                double curr_mag_of_bfield = rmag(Bx(BFx,i,y,j),By(BFy,i,y,j),Bz(BFz,i,y,j));
                fprintf(bfield_test_file, "%d %d %lf\n", i, j, curr_mag_of_bfield);               
            }
        }
    }
    else if(z_plane == 1){
        printf("Creating data file for a heat map for an x-y plane at z = %lf!\n", z);
        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                double curr_mag_of_bfield = rmag(Bx(BFx,i,j,z),By(BFy,i,j,z),Bz(BFz,i,j,z));               
                fprintf(bfield_test_file, "%d %d %lf\n", i, j, curr_mag_of_bfield);            
            }
        }
    }
    else{
        printf("No plane selected! Exiting program.\n");
        exit(0);
    }

    fclose(bfield_test_file);
    printf("Data file has been created! Exiting c code!\n");

}

int main(void){
    const char filename[] = "Bfield.dat";
    magnetic_field_test(filename);
    return 0;
}
