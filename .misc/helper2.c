int fun1(int x) { printf("Fun1\n"); return x; } 
int fun2(int x) { printf("Fun2\n"); return 2*x; } 



// #include <stdio.h> 
// // A normal function with an int parameter 
// // and void return type 
// int eloss1(int a) 
// { 
//     return a;
// } 
// int eloss2(int a){

//   return 2*a;  
// }
// int (*fun_ptr)(int);
  
// int main() 
// { 
//     int (*fun_ptr)(int) = &eloss1;
//     int (*fun_ptr)(int) = &eloss2;
//     int part = 2;
//     if(part == 1){
//     // fun_ptr is a pointer to function fun()  
    
//     int (*fun_ptr)(int) = &eloss1; }
//     else{
//     printf("Got here\n");
//     int (*fun_ptr)(int) = &eloss2;  
//     }
  
//     /* The above line is equivalent of following two 
//        void (*fun_ptr)(int); 
//        fun_ptr = &fun;  
//     */
  
//     // Invoking fun() using fun_ptr  
//     printf("%d\n",(*fun_ptr)(15));
//     return 0; 
// } 