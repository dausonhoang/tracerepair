#include <stdio.h>
#include <stdlib.h>
#include <string.h>		// for memset, memcmp
#include "erasure_code.h"
#include "types.h"

void gf_vect_mul_init(unsigned char c, unsigned char *tbl)
{
	unsigned char c2 = (c << 1) ^ ((c & 0x80) ? 0x1d : 0);	//Mult by GF{2}
	unsigned char c4 = (c2 << 1) ^ ((c2 & 0x80) ? 0x1d : 0);	//Mult by GF{2}
	unsigned char c8 = (c4 << 1) ^ ((c4 & 0x80) ? 0x1d : 0);	//Mult by GF{2}
	unsigned char c3, c5, c6, c7, c9, c10, c11, c12, c13, c14, c15;
	unsigned char c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, c28, c29, c30,
	    c31;

	c3 = c2 ^ c;
	c5 = c4 ^ c;
	c6 = c4 ^ c2;
	c7 = c4 ^ c3;

	c9 = c8 ^ c;
	c10 = c8 ^ c2;
	c11 = c8 ^ c3;
	c12 = c8 ^ c4; 
	c13 = c8 ^ c5;
	c14 = c8 ^ c6;
	c15 = c8 ^ c7;

	tbl[0] = 0;
	tbl[1] = c;
	tbl[2] = c2;
	tbl[3] = c3;
	tbl[4] = c4;
	tbl[5] = c5;
	tbl[6] = c6;
	tbl[7] = c7;
	tbl[8] = c8;
	tbl[9] = c9;
	tbl[10] = c10;
	tbl[11] = c11;
	tbl[12] = c12;
	tbl[13] = c13;
	tbl[14] = c14;
	tbl[15] = c15;

	c17 = (c8 << 1) ^ ((c8 & 0x80) ? 0x1d : 0);
	 /*if c8 > 128 -> XOR 29. 
    Else XOR 0

    case 1: c8 << 1 XOR 29 
    case 2: c8 << 1 XOR 0*/ 
	c18 = (c17 << 1) ^ ((c17 & 0x80) ? 0x1d : 0);	//Mult by GF{2}//29 = 0001 1101
	c19 = c18 ^ c17;
	c20 = (c18 << 1) ^ ((c18 & 0x80) ? 0x1d : 0);	//Mult by GF{2}
	c21 = c20 ^ c17;
	c22 = c20 ^ c18;
	c23 = c20 ^ c19;
	c24 = (c20 << 1) ^ ((c20 & 0x80) ? 0x1d : 0);	//Mult by GF{2}
	c25 = c24 ^ c17;
	c26 = c24 ^ c18;
	c27 = c24 ^ c19;
	c28 = c24 ^ c20;
	c29 = c24 ^ c21;
	c30 = c24 ^ c22;
	c31 = c24 ^ c23;

	tbl[16] = 0;
	tbl[17] = c17;
	tbl[18] = c18;
	tbl[19] = c19;
	tbl[20] = c20;
	tbl[21] = c21;
	tbl[22] = c22;
	tbl[23] = c23;
	tbl[24] = c24;
	tbl[25] = c25;
	tbl[26] = c26;
	tbl[27] = c27;
	tbl[28] = c28;
	tbl[29] = c29;
	tbl[30] = c30;
	tbl[31] = c31;

    for (int i = 0; i < 32; i++){
        printf("%u ", tbl[i]);
    }
}

int main(){

    unsigned char *g_tbls;
    //gf_vect_mul_init(162, g_tbls);
  
    /*unsigned char c;
    c = 162;
    unsigned char c2;
    c2 = (c << 1) ^ ((c & 0x80) ? 0x1d : 0);	//Mult by GF{2}
	unsigned char c4;
    c4 = (c2 << 1) ^ ((c2 & 0x80) ? 0x1d : 0);	//Mult by GF{2}
	unsigned char c8;
    c8 = (c4 << 1) ^ ((c4 & 0x80) ? 0x1d : 0);	//Mult by GF{2}

    

    printf("c2: %u", c2);*/
    return 0;
}