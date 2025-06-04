  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <time.h>
  #include <string.h>
  #include <assert.h>
  #include <stdbool.h>
  #include <float.h>


int main(int argc, char ** argv) {
    int current_run_number = 100000000;
    FILE *ptfileout; // output file
    int size = sizeof(current_run_number);
    char *fileout_name, number_str[size];
    sprintf(number_str,"%d", current_run_number);
    printf("%s\n",number_str);
    fileout_name = strcat(number_str,"_output.dat");
    printf("%s\n",fileout_name);
    ptfileout = fopen(fileout_name, "w");
    printf("Opened file successfully\n");
    return 0;
}