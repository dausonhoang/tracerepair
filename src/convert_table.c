#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"repair_scheme_tbl.h"

/*void decToHexa(unsigned char n)
{  
    char hexaDeciNum[100];
    

    if(n == 0){
        printf("0x00, ");
    }

    int i = 0;
    while (n != 0) {
        int temp = 0;
        temp = n % 16;
        if (temp < 10) {
            hexaDeciNum[i] = temp + 48;
            i++;
        }
        else {
            hexaDeciNum[i] = temp + 55;
            i++;
        }
        n = n / 16;
    }
    
    printf("0x");
    for (int j = i - 1; j >= 0; j--){
        
        printf("%c", hexaDeciNum[j]);
    }
    printf(", ");
}
*/

unsigned char *binRep(unsigned char n, unsigned char t){

    unsigned char* vec = (unsigned char*)malloc(8*sizeof(unsigned char));
	if (vec == NULL){
		printf("Error: Allocate fail vec-binrep");
		free(vec);
	}

    for (int i = 1; i <= t; i++){
        vec[t-i] = n & 1;
        n = n >> 1;
    }
    return vec;
}


int main(){
    int i, j, z;
    unsigned char *p;
    int vi;


    for (i = 0; i < 9; i++){ 
        printf("{");
        for (j = 0; j < 9; j++){
            printf("{ ");
            for(z = 0; z < 9; z++){
                p = binRep(h_rtbl[i][j][z], h_htbl[i][j][0]);
                for( int l = 0; l < h_htbl[i][j][0]; l++){
                    printf("%u, ", p[l]);
                }
                printf("| ");
            }
            printf("}, \n");
        }
        printf(" },\n");
    }
    printf("\n");

    
    return 0;

}

