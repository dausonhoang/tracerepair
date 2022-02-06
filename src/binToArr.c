#include<stdio.h>
#include<stdlib.h>

int main(){
    unsigned char sum;
    unsigned char *binArr;
    int rem; 
    int i, j;
    i = 0;

    binArr=malloc(8);

    scanf("%d", &sum);
    while(sum > 0){
        rem = sum % 10;
        binArr[i] = rem;
        sum = sum/10;
        i++;
    }
    for(j=i-1; j>=0; j--){
        printf("%d-", binArr[j]);
    }
    free(binArr);
    return 0;
}