#include "ec_base.h"
#include <string.h>	
#include <stdio.h>
#include <stdlib.h> 

unsigned char gf_inv(unsigned char a)
{
#ifndef GF_LARGE_TABLES
	if (a == 0)
		return 0;
	return gff_base[255 - gflog_base[a & 0xff] & 0xff];
	// if a = 1111.1111 & 1111.1111 => 1111.1111 
#else
	return gf_inv_table_base[a];
#endif
}	

int main(int argc, char *argv[])
{
	
    int m, k,i,j;

    m = 9;
	k = 6;
	for (i = k; i < m; i++)
	{
		for (j = 0; j < k; j++) 
		
			printf("%u ",gf_inv(i ^ j));
		printf("\n");
	
	}
	
    return 0;
}