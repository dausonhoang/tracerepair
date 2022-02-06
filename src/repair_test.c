#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
unsigned char *binGap(int t, int m){

    int q; 
    int i, log; 
    static unsigned char bin[8];

    if (m < 0 || m > pow(2, (t))-1){
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
   unsigned char buffs[9];
	FILE *myFile;
    myFile = fopen("RandomCodeword96.txt", "r");

    if (myFile == NULL){
        printf("Error Reading File\n");
        exit (0);
    }

    for (int a = 0; a < 9; a++)
    {
        fscanf(myFile, "%d,", &buffs[a]);
    }
     for (int a = 0; a < 9; a++)
    {
        printf("%u", buffs[a]);
    }
	fclose(myFile);


    return 0;
}