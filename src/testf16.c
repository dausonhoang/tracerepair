#include "ec_base.h"
#include <string.h>	
#include <stdio.h>
#include <stdlib.h> 


unsigned char gf_mul(unsigned char a, unsigned char b)
{
#ifndef GF_LARGE_TABLES
	int i;

	if ((a == 0) || (b == 0))
		return 0;

	return gff_base_f16[(i = gflog_base_f16[a] + gflog_base_f16[b]) > 254 ? i - 255 : i];
#else
	return gf_mul_table_base[b * 256 + a];
#endif

}

int main()
{
    int a, b;

   printf("enter a: ");
   scanf("%u", &a);

   printf("enter b: ");
   scanf("%u", &b);

    printf("Output %u", gf_mul(a, b));
    return 0;

}