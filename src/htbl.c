#include<stdlib.h>
#include<stdio.h>
#include"repair_scheme_tbl.h"

int main(){
    int (*htbl)[9][9][9] = &h_htblb;
    for(int a=0; a<9; a++){
		for(int b=0; b<9; b++){
			for(int c = 0; c<9; c++){
				printf("%d \n", htbl[a][b][c]);
			}
		}
	}
    return 0;
}