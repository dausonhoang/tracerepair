#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "erasure_code.h"
#include "ec_base.h"
#include<time.h>

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

void gf_vect_mul_init(unsigned char c, unsigned char *tbl)
{
	unsigned char c2 = (c << 1) ^ ((c & 0x80) ? 0x1d : 0);	//Mult by GF{2}
	unsigned char c4 = (c2 << 1) ^ ((c2 & 0x80) ? 0x1d : 0);	//Mult by GF{2}
	unsigned char c8 = (c4 << 1) ^ ((c4 & 0x80) ? 0x1d : 0);	//Mult by GF{2}

#if __WORDSIZE == 64 || _WIN64 || __x86_64__
	unsigned long long v1, v2, v4, v8, *t;
	unsigned long long v10, v20, v40, v80;
	unsigned char c17, c18, c20, c24;

	t = (unsigned long long *)tbl;

	v1 = c * 0x0100010001000100ull;
	v2 = c2 * 0x0101000001010000ull;
	v4 = c4 * 0x0101010100000000ull;
	v8 = c8 * 0x0101010101010101ull;

	v4 = v1 ^ v2 ^ v4;
	t[0] = v4;
	t[1] = v8 ^ v4;

	c17 = (c8 << 1) ^ ((c8 & 0x80) ? 0x1d : 0);	//Mult by GF{2}
	c18 = (c17 << 1) ^ ((c17 & 0x80) ? 0x1d : 0);	//Mult by GF{2}
	c20 = (c18 << 1) ^ ((c18 & 0x80) ? 0x1d : 0);	//Mult by GF{2}
	c24 = (c20 << 1) ^ ((c20 & 0x80) ? 0x1d : 0);	//Mult by GF{2}

	v10 = c17 * 0x0100010001000100ull;
	v20 = c18 * 0x0101000001010000ull;
	v40 = c20 * 0x0101010100000000ull;
	v80 = c24 * 0x0101010101010101ull;

	v40 = v10 ^ v20 ^ v40;
	t[2] = v40;
	t[3] = v80 ^ v40;

#else // 32-bit or other 16 bit
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

	c17 = (c8 << 1) ^ ((c8 & 0x80) ? 0x1d : 0);	//Mult by GF{2}
	c18 = (c17 << 1) ^ ((c17 & 0x80) ? 0x1d : 0);	//Mult by GF{2}
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

#endif //__WORDSIZE == 64 || _WIN64 || __x86_64__
}

void gf_gen_rs_matrix(unsigned char *a, int m, int k)
{
	int i, j;
	unsigned char p, gen = 1;

	memset(a, 0, k * m);
	for (i = 0; i < k; i++)
		a[k * i + i] = 1;

	for (i = k; i < m; i++) {
		p = 1;
		for (j = 0; j < k; j++) {
			a[k * i + j] = p;
			p = gf_mul(p, gen);
		}
		gen = gf_mul(gen, 2);
	}
}

void ec_init_tables(int k, int rows, unsigned char *a, unsigned char *g_tbls)
{
	int i, j;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < k; j++) {
			gf_vect_mul_init(*a++, g_tbls); // gf_vect_mul_inti: Calculates const table gftbl in GF(2^8) from single input A
			g_tbls += 32;
		}
    }
}

void ec_encode_data_base(int len, int srcs, int dests, unsigned char *v,
			unsigned char **src, unsigned char **dest)
			
//Len: Length of each block of data (vector) of source or dest data.
//srcs (k): The number of vector sources or rows in the generator matrix for coding
//dests (rows): The number of output vectors to concurrently encode/decode.
//v(g_tbls): Pointer to array of input tables generated from coding coefficients in ec_init_tables(). 
//Must be of size 32*k*rows
//src (data): Array of pointers to source input buffers.
//dest (coding): Array of pointers to coded output buffers.

{
	int i, j, l;
	unsigned char s;

	for (l = 0; l < dests; l++) {
		for (i = 0; i < len; i++) {
			s = 0;
			for (j = 0; j < srcs; j++)
				s ^= gf_mul(src[j][i], v[j * 32 + l * srcs * 32 + 1]);
				// v is pointer to g_tbls where they pre-compute 32 elements
				// g_tbls is where they generate each table that has 32 elements. If we have n, k > 32 they keep computing other 32 elements table 

			dest[l][i] = s;
		}
	}
}

// Generate Random errors
/*static void gen_err_list(unsigned char *src_err_list,
			 unsigned char *src_in_err, int *pnerrs, int *pnsrcerrs, int k, int m)
{
	int i, err;
	int nerrs = 0, nsrcerrs = 0;

	for (i = 0, nerrs = 0, nsrcerrs = 0; i < m && nerrs < m - k; i++) {
		err = 1 & rand();
		src_in_err[i] = err;
		if (err) {
			src_err_list[nerrs++] = i;
			if (i < k) {
				nsrcerrs++;
			}
		}
	}
	if (nerrs == 0) {	// should have at least one error
		while ((err = (rand() % KMAX)) >= m) ;
		src_err_list[nerrs++] = err;
		src_in_err[err] = 1;
		if (err < k)
			nsrcerrs = 1;
	}
	*pnerrs = nerrs;
	*pnsrcerrs = nsrcerrs;
	return;
}*/


#define TEST_LEN 20
#define TEST_SOURCES 127
#define MAX TEST_SOURCES

int main(int argc, char *argv[])
{

    int m = 14; 
    int k = 9;


    int rows = m - k;
    unsigned char *encode_matrix;
    unsigned char *g_tbls;
	unsigned char *buffs[TEST_SOURCES], *temp_buffs[TEST_SOURCES];
	void *buf;
    int i, j;

	//allocate array
	for (i = 0; i < TEST_SOURCES; i++){
		buffs[i] = buf;
	}

	for (i = 0; i < TEST_SOURCES; i++){
		temp_buffs[i] = buf;
	}

	
    encode_matrix = malloc(MAX * MAX);
	g_tbls = malloc(MAX * MAX * 32);
	if (encode_matrix == NULL || g_tbls == NULL) {
		printf("Test failure! Error with malloc\n");
		return -1;
	}

	/*
	//Data rand
	unsigned char a[6][10] = {
		{12, 34, 45, 74, 2, 48, 34, 5, 75, 10}, 
		{76, 87, 8, 23, 55, 19, 98, 39, 85, 4},
		{86, 38, 41, 36, 30, 14, 47, 25, 1, 9},
		{61, 85, 28, 99, 52, 77, 44, 7, 9, 23},
		{75, 24, 52, 65, 97, 2, 13, 54, 65, 91},
		{22, 76, 3, 94, 25, 65, 24, 54, 68, 11}
	};

	for (i = 0; i < TEST_SOURCES; i++)
	{
		for (j = 0; j < sizeof(a); j++)
		{
			buffs[i] = (unsigned char*)a[j];
		}
		break;
	}

	for (i = 0; i < sizeof(buffs); i++)
	{
		printf("%u ", *buffs[i]);
	}
	*/
	
	for (i = 0; i < k; i++)
	{ 
		for (j = 0; j < TEST_LEN; j++)
		{
			buffs[i][j] = rand();
			printf("%u ", *buffs[j]);

			if(j == 30)
			printf("\n");
		}
	}

	
    gf_gen_rs_matrix(encode_matrix, m, k);
	printf("The G matrix is: \n");

	for(i = 0 ; i < k*m ; i++)
	{
		printf("%u ",encode_matrix[i]);
		if((i % k) == k-1)
		printf("\n");
	}
	
    ec_init_tables(k, m - k, &encode_matrix[k * k], g_tbls);
	printf("The initial table: \n");
	for (i = 0; i < rows*k; i++)
	{
		printf("%u ", g_tbls[i]);
		if (i == 32)
		printf("\n");
	}
	printf("\n");

	ec_encode_data_base(TEST_LEN, k, m - k, g_tbls, buffs, &buffs[k]);

	for(i = 0; i < k*TEST_LEN ; i++)
	{
		printf("%u ", *buffs[i]);
	}


    free(encode_matrix);
    free(g_tbls);
    return 0;


}
