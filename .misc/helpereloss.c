// A simple C program to show function pointers as parameter 
#include <stdio.h>
#include "helper2.h"  
// Two simple functions
  
// A function that receives a simple function 
// as parameter and calls the function 
//void wrapper( int (*fun)(int)) 
//{ 
//    printf("%d\n",(fun(1))); 
//} 
  
int main() 
{ 
  //  wrapper(fun1); 
//    wrapper(fun2); 
    int (*help)(int);
    int part = 1;

    if(part == 1){
    help = fun1;}
    else{
    help = fun2;
    }
    printf("%d\n",help(3));


    return 0; 
}



















// #include <stdio.h> 
// #include <stdlib.h>
// // A normal function with an int parameter 
// // and void return type 
// int fun(int a) 
// { 
//     //printf("Value of a is %d\n", a); 
//     return a;
// } 
  
// int fun2(int a) 
// { 
//     //printf("Value of a is %d\n", a); 
//     return 2*a;
// } 

// int (*fun_ptr)(int);
  
// int main() 
// {  
//     int part = 1;
//     if(part == 1){
//     int (*fun_ptr)(int) = fun2;  // & removed 
//   }
//     printf("%d\n",fun_ptr(10));  // * removed 
  
//     return 0; 
// }