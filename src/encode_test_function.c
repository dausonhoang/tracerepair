#include <limits.h>
#include <string.h>	
#include <stdio.h>
#include <stdlib.h>
#include "erasure_code.h"
#include "ec_base.h"		// for GF tables


unsigned char gf_mul(unsigned char a, unsigned char b)
{
	#ifndef GF_LARGE_TABLES
	int i;

	if ((a == 0) || (b == 0))
		return 0;

	return gff_base[(i = gflog_base[a] + gflog_base[b]) > 254 ? i - 255 : i];
#else
	return gf_mul_table_base[b * 256 + a];
#endif
}

int main(){
	unsigned char c;
	c = gf_mul(74, 1);
	printf("%u ", c);
	return 0;



}	
