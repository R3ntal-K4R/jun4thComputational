#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#define clear "\033[2K"
#define startline "\033[<10>D"
#define restore "\033[u"
#define move "\033[<0>;<0>H"

int main(){
	while(1){
    for(int i=0; i<10; i++){
        printf("=");
        fflush(stdout);
        usleep(100000);
}
printf(">");
printf(clear);
printf("\r");}
printf("\n");
    
}
