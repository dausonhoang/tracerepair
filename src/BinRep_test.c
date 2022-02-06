#include <stdio.h>
#include <stdlib.h>

int* binRep(unsigned int n, int t){
   
    int* vec = (int*) malloc(t);
   
    for (int i = 1; i <= t; i++){
        vec[t-i] = n & 1;
        n = n >> 1;
    }
    return vec;
}

int main()
{
    int t = 8;
    int n = 7; 

    int* rep = binRep(n, t);
    for (int i = 0; i < 8; i++) printf(" %d ", rep[i]);
   
    return 0;
}