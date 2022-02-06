#include <stdlib.h>
#include<stdio.h>
#include<math.h>


int binRep(unsigned char N) 
{ 
    // To store the binary number 
    int B_Number = 0; 
    int c;
    int count = 0; 

    while (N != 0) { 
        int rem = N % 2; 
        c = pow(10, count); 

        B_Number += rem * c; 

        N /= 2; 
        count++; 
    }
    return B_Number;
}

struct{
	unsigned int traceSum:1;

}bit_field;

unsigned int binSum(int x, int y){

    unsigned char sum;
    unsigned char binArr[8];
    int rem, sum1; 
    int i, a, b;
    bit_field.traceSum = 0;
    
    sum = x^y;
    printf("\na^b: %u", binRep(sum));
    sum1 = binRep(sum);

    a = 0;
    while(sum1 > 0){
        rem = sum1 % 10;
        binArr[a] = rem;
        sum1 = sum1/10;
        a++;
    }

    printf("\nbin Arr: ");
    for(b=a-1; b>=0; b--){
        printf(" %u ", binArr[b]);
    }
    
	for(i = 0; i < 8 ; i++){
		bit_field.traceSum ^= binArr[i];
        if(binArr[i]!=0 || binArr[i]!=1)
            break;
	}
   
    return bit_field.traceSum;
}

unsigned char *binGap(int t, int m){
    int q; 
    int i, log; 
    static unsigned char bin[8];

    if (m < 0 || m > pow(2, (t-1))){
        printf("Number not in range [0..(q^t-1)]");
    }

    for (i = 0; i < (t); i++){
        bin[i] = 0;
    }
    while ( m > 1){
        log = log2(m);
        bin[t-log-1] = 1;
        m = m - pow(2, log);
    }
    
    if( m == 1){
        bin[t-1] = 1;
    }
    return bin;

}


int main(){

    int a = 200;
    int b = 153;
    int t, m;

    t = 5;
    m = 15;

    unsigned char *p; 
    
    p = binGap(t, m);
    
    bit_field.traceSum = binSum(a,b);
    printf("\nSum : %u", bit_field.traceSum);
    printf("\nSize of bit_field: %d\n", sizeof(bit_field));
    
    printf("binGAP:");
    
    
    //for(int i=t-1; i>=0; i--){
    for (int i = 0; i < t; i++){
        printf("%u ", p[i]);
    }
    
    return 0;

}