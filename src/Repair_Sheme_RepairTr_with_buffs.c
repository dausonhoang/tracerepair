/**********************************************************************
  Copyright(c) 2011-2015 Intel Corporation All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
    * Neither the name of Intel Corporation nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		// for memset, memcmp
#include "erasure_code.h"
#include "types.h"
#include<time.h>
#include<math.h>
#include "repair_scheme_tbl_9_6.h"
#include <assert.h>

//add

#define TEST_SOURCES 127
#define TEST_LEN 100000
#define BIN_LEN 8
//100.000 = 1.2->1.5s
//0.015
//end add

#define TEST_SIZE (TEST_LEN/2)


#ifndef TEST_SOURCES
# define TEST_SOURCES  127
#endif
#ifndef RANDOMS
# define RANDOMS 100
#endif

#define MMAX TEST_SOURCES
#define KMAX TEST_SOURCES

#define EFENCE_TEST_MIN_SIZE 16
#define EFENCE_TEST_MAX_SIZE EFENCE_TEST_MIN_SIZE + 0x100

#ifdef EC_ALIGNED_ADDR
// Define power of 2 range to check ptr, len alignment
# define PTR_ALIGN_CHK_B 0
# define LEN_ALIGN_CHK_B 0	// 0 for aligned only
#else
// Define power of 2 range to check ptr, len alignment
# define PTR_ALIGN_CHK_B 32
# define LEN_ALIGN_CHK_B 32	// 0 for aligned only
#endif

#ifndef TEST_SEED
#define TEST_SEED 11
#endif

#include <limits.h>
#include <string.h>		// for memset
#include "erasure_code.h"
#include "ec_base.h"		// for GF tables'#include<sys/time.



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



unsigned char gf_inv(unsigned char a)
{
#ifndef GF_LARGE_TABLES
	if (a == 0)
		return 0;

	return gff_base[255 - gflog_base[a]];
#else
	return gf_inv_table_base[a];
#endif
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

void gf_gen_cauchy1_matrix(unsigned char *a, int m, int k)
{
	int i, j;
	unsigned char *p;

	// Identity matrix in high position
	memset(a, 0, k * m);
	for (i = 0; i < k; i++)
		a[k * i + i] = 1;

	// For the rest choose 1/(i + j) | i != j
	p = &a[k * k];
	for (i = k; i < m; i++)
		for (j = 0; j < k; j++)
			*p++ = gf_inv(i ^ j);

}

int gf_invert_matrix(unsigned char *in_mat, unsigned char *out_mat, const int n)
{
	int i, j, k;
	unsigned char temp;

	// Set out_mat[] to the identity matrix
	for (i = 0; i < n * n; i++)	// memset(out_mat, 0, n*n)
		out_mat[i] = 0;

	for (i = 0; i < n; i++)
		out_mat[i * n + i] = 1;

	// Inverse
	for (i = 0; i < n; i++) {
		// Check for 0 in pivot element
		if (in_mat[i * n + i] == 0) {
			// Find a row with non-zero in current column and swap
			for (j = i + 1; j < n; j++)
				if (in_mat[j * n + i])
					break;

			if (j == n)	// Couldn't find means it's singular
				return -1;

			for (k = 0; k < n; k++) {	// Swap rows i,j
				temp = in_mat[i * n + k];
				in_mat[i * n + k] = in_mat[j * n + k];
				in_mat[j * n + k] = temp;

				temp = out_mat[i * n + k];
				out_mat[i * n + k] = out_mat[j * n + k];
				out_mat[j * n + k] = temp;
			}
		}

		temp = gf_inv(in_mat[i * n + i]);	// 1/pivot
		for (j = 0; j < n; j++) {	// Scale row i by 1/pivot
			in_mat[i * n + j] = gf_mul(in_mat[i * n + j], temp);
			out_mat[i * n + j] = gf_mul(out_mat[i * n + j], temp);
		}

		for (j = 0; j < n; j++) {
			if (j == i)
				continue;

			temp = in_mat[j * n + i];
			for (k = 0; k < n; k++) {
				out_mat[j * n + k] ^= gf_mul(temp, out_mat[i * n + k]);
				in_mat[j * n + k] ^= gf_mul(temp, in_mat[i * n + k]);
			}
		}
	}
	return 0;
}

//
// gftbl(A) = {A{00}, A{01}, A{02}, ... , A{0f} 15 }, {A{00}, A{10}, A{20}, ... , A{f0} 240 }

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

#endif //__WORDSIZE == 64 || _WIN64 || __x86_64__
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

	for (l = 0; l < dests; l++) { // l = 1
		for (i = 0; i < len; i++) { // i = 0
			s = 0;
			for (j = 0; j < srcs; j++){
				s^= gf_mul(src[j][i], v[j * 32 + l * srcs * 32 + 1]);
			}
			dest[l][i] = s;
			//printf("s: %u \n", s);
		}
	}
}

void ec_encode_data_update_base(int len, int k, int rows, int vec_i, unsigned char *v,
				unsigned char *data, unsigned char **dest)
{
	int i, l;
	unsigned char s;

	for (l = 0; l < rows; l++) {
		for (i = 0; i < len; i++) {
			s = dest[l][i];
			s ^= gf_mul(data[i], v[vec_i * 32 + l * k * 32 + 1]);

			dest[l][i] = s;
		}
	} //&buff[i]
}


struct slver {
	unsigned short snum;
	unsigned char ver;
	unsigned char core;
};

// Version info
struct slver gf_vect_mul_init_slver_00020035;
struct slver gf_vect_mul_init_slver = { 0x0035, 0x02, 0x00 };

struct slver ec_encode_data_base_slver_00010135;
struct slver ec_encode_data_base_slver = { 0x0135, 0x01, 0x00 };

struct slver gf_vect_mul_base_slver_00010136;
struct slver gf_vect_mul_base_slver = { 0x0136, 0x01, 0x00 };

struct slver gf_vect_dot_prod_base_slver_00010137;
struct slver gf_vect_dot_prod_base_slver = { 0x0137, 0x01, 0x00 };

struct slver gf_mul_slver_00000214;
struct slver gf_mul_slver = { 0x0214, 0x00, 0x00 };

struct slver gf_invert_matrix_slver_00000215;
struct slver gf_invert_matrix_slver = { 0x0215, 0x00, 0x00 };

struct slver gf_gen_rs_matrix_slver_00000216;
struct slver gf_gen_rs_matrix_slver = { 0x0216, 0x00, 0x00 };

struct slver gf_gen_cauchy1_matrix_slver_00000217;
struct slver gf_gen_cauchy1_matrix_slver = { 0x0217, 0x00, 0x00 };

void ec_encode_data(int len, int srcs, int dests, unsigned char *v,
		    unsigned char **src, unsigned char **dest)
{
	ec_encode_data_base(len, srcs, dests, v, src, dest);
}


typedef unsigned char u8;

void dump(unsigned char *buf, int len)
{
	int i;
	for (i = 0; i < len;) {
		printf(" %u", 0xff & buf[i++]);
		if (i % 32 == 0)
		printf("\n");
	}
	printf("\n");
}

void dump_matrix(unsigned char **s, int k, int m)
{
	int i, j;
	for (i = 0; i < k; i++) {
		for (j = 0; j < m; j++) {
			printf(" %u", s[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void dump_u8xu8(unsigned char *s, int k, int m)
{
	int i, j;
	for (i = 0; i < k; i++) {
		for (j = 0; j < m; j++) {
			printf(" %u", 0xff & s[j + (i * m)]);
		}
		printf("\n");
	}
	printf("\n");
}

// Generate Random errors
static void gen_err_list(unsigned char *src_err_list,
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
	*pnsrcerrs = nsrcerrs;// number of data error
	return;
}

#define NO_INVERT_MATRIX -2
// Generate decode matrix from encode matrix
static int gf_gen_decode_matrix(unsigned char *encode_matrix,
				unsigned char *decode_matrix,
				unsigned char *invert_matrix,
				unsigned int *decode_index,
				unsigned char *src_err_list,
				unsigned char *src_in_err,
				int nerrs, int nsrcerrs, int k, int m)
{
	int i, j, p;
	int r;
	unsigned char *backup, *b, s;
	int incr = 0;

	b = malloc(MMAX * KMAX);
	//MALLOC: or �memory allocation� method in C is used to dynamically allocate
	//a single large block of memory with the specified size.
	//It returns a pointer of type void which can be cast into a pointer of any form.
	// It initializes each block with default garbage value

	backup = malloc(MMAX * KMAX);

	if (b == NULL || backup == NULL) {
		printf("Test failure! Error with malloc\n");
		free(b);
		free(backup);
		//is used to dynamically de-allocate the memory.
		//The memory allocated using functions malloc() is not de-allocated on their own.
		//Hence the free() method is used. helps to reduce wastage of memory by freeing it
		return -1;
	}

	// Construct matrix b by removing error rows
	//printf("B AND BACKUP MATRIX: \n");
	for (i = 0, r = 0; i < k; i++, r++) {
		while (src_in_err[r])
			r++;
		for (j = 0; j < k; j++) {
			b[k * i + j] = encode_matrix[k * r + j];
			backup[k * i + j] = encode_matrix[k * r + j];
			/*printf("%u ", b[k * i + j]);
			if ((j%k)==k-1)
			printf("\n");*/
		}
		decode_index[i] = r; //non erased data location
	}
	incr = 0;

	while (gf_invert_matrix(b, invert_matrix, k) < 0) {
		if (nerrs == (m - k)) {
			free(b);
			free(backup);
			printf("BAD MATRIX\n");
			return NO_INVERT_MATRIX;
		}
		incr++;
		memcpy(b, backup, MMAX * KMAX);
		for (i = nsrcerrs; i < nerrs - nsrcerrs; i++) {
			if (src_err_list[i] == (decode_index[k - 1] + incr)) {
				// skip the erased parity line
				incr++;
				continue;
			}
		}
		if (decode_index[k - 1] + incr >= m) {
			free(b);
			free(backup);
			printf("BAD MATRIX\n");
			return NO_INVERT_MATRIX;
		}
		decode_index[k - 1] += incr;
		for (j = 0; j < k; j++)
			b[k * (k - 1) + j] = encode_matrix[k * decode_index[k - 1] + j];

	};

	//non parity decoding. eg recover data that has index < k
	for (i = 0; i < nsrcerrs; i++) {
		for (j = 0; j < k; j++) {
			decode_matrix[k * i + j] = invert_matrix[k * src_err_list[i] + j];
		}
	}
	/* src_err_list from encode_matrix * invert of b for parity decoding */
	for (p = nsrcerrs; p < nerrs; p++) { //nserrs: number of erased data unit , nerrs: Erased index
		for (i = 0; i < k; i++) {
			s = 0;
			for (j = 0; j < k; j++)
				s ^= gf_mul(invert_matrix[j * k + i],
					    encode_matrix[k * src_err_list[p] + j]);

			decode_matrix[k * p + i] = s;
		}// b is parity check matrix
	}
	free(b);
	free(backup);
	return 0;
}
/********************************************************************************/

//version with RepairTr *buffs
//This version of repair_trace_optimised use 1-D array approach and also optimise the data layout to leverage caching
//repair_trace / repair_trace_optimised 3D / repair_trace_optimised 1D : 42x / 35x / 34x (for TEST_LEN = 100K, for 1M, it fluctuates)
//1D approach can run with larger TEST_LEN, etc. 1,000,000, while 3D cannot

unsigned char* repair_trace_optimised(int block, int j, unsigned char **buffs){
	int a, b, c, test_codeword; // for FOR loop
  	int i, s, vi;
	unsigned char bi[block];
	unsigned char *vec;
	unsigned char result;

	//time_t st = clock();



	unsigned char* RepairTr = (unsigned char*) malloc(sizeof(unsigned char)*block*TEST_LEN);
	unsigned char* TargetTr = (unsigned char*) malloc(sizeof(unsigned char)*8*TEST_LEN);
	unsigned char* rev = (unsigned char*) malloc(sizeof(unsigned char)*TEST_LEN);


	if (!RepairTr){
		printf("Cannot allocate memory for RepairTr"); exit(1);
	}
		if (!TargetTr){
		printf("Cannot allocate memory for TargetTr"); exit(1);
	}

	/*
	time_t fi = clock();
	double malloc_time = ((double)(fi - st) / CLOCKS_PER_SEC) /1.0;
    printf("Malloc time: %f seconds\n", malloc_time);
	*/

	if (!rev){
		printf("Cannot allocate memory for rev"); exit(1);
	}

    //Calculate bandwidth using H table
	//printf("Bandwidth: ");
    for (i = 0; i < block; i++){
			if (i != j){
				bi[i] = h_htbl[i][j][0];
			}
	}

		//Calculate traces to send to Node j (H table)

		//unsigned char *buffs_p;
		unsigned char buffs_p, u;
		unsigned char table[72];

		//printf("\nRepair Traces: \n");


		for (i = 0; i < block ; i++){
			if( i != j){
				//use 1-D array table to reduce the access to h_htbl
				for (u = 0; u < 9; u++){
					table[u] = h_htbl[i][j][u];
				}
                for (test_codeword = 0; test_codeword < TEST_LEN; test_codeword++){
                    buffs_p = buffs[i][test_codeword];
                    for (a = 0; a < bi[i]; a++){
						RepairTr[i*TEST_LEN+test_codeword] |= (parity[table[a+1] & buffs_p]) << a;
					}
				}
			}
		}

	//Recover node generates Column traces using R table
	//printf("Column Tr: \n");
    for (i = 0; i < block; i++){
		if (i != j){
			for (u = 0; u < 8; u++){
					table[u] = h_rtbl[i][j][u];
			}
			//p = RepairTr[i];
			//p = &(RepairTr[i*block*TEST_LEN]);
			for (test_codeword = 0; test_codeword < TEST_LEN; test_codeword++){
        		for (s = 0; s < 8; s++){
					TargetTr[test_codeword * 8 + s] ^= (parity[table[s] & RepairTr[i*TEST_LEN+test_codeword]]);
				}
    		}
		}
	}

	//printf("\n\n");
	//Assume that we obtained D table assigned as dtbl
	for (u = 0; u < 8; u++){
		table[u] = h_dtbl[j][u];
	}
	for (test_codeword = 0; test_codeword < TEST_LEN; test_codeword++){
		rev[test_codeword] = 0; //0*Z(q)
		for(s = 0; s < 8; s++){
			//rev[test_codeword] ^= TargetTr[s][test_codeword]*table[s];
			//rev[test_codeword] ^= TargetTr[s*TEST_LEN+test_codeword]*table[s];
			rev[test_codeword] ^= TargetTr[test_codeword*8+s]*table[s];			//better data layout
		}
	}
	return rev;
}

void precompute() {
    int i;
    parity[0] = 0;
    for (i = 1; i < 256; i++) {
        parity[i] = parity[i >> 1] ^ (i & 1);
    }
}

int main(int argc, char *argv[])
{
	int re = 0;
	int i, j, p, rtest, m, k;
	int nerrs, nsrcerrs;
	void *buf;
	unsigned int decode_index[MMAX];
	unsigned char *temp_buffs[TEST_SOURCES], *buffs[TEST_SOURCES];
	unsigned char *encode_matrix, *decode_matrix, *invert_matrix, *g_tbls;
	unsigned char src_in_err[TEST_SOURCES], src_err_list[TEST_SOURCES];
	unsigned char *recov[TEST_SOURCES];

	int rows, align, size;
	unsigned char *efence_buffs[TEST_SOURCES];
	unsigned int offset;
	u8 *ubuffs[TEST_SOURCES];
	u8 *temp_ubuffs[TEST_SOURCES];

	printf("erasure_code_test: %dx%d \n", TEST_SOURCES, TEST_LEN);
	//srand(TEST_SEED);
	srand(time(0));		//to make random data change with different calls


	// Allocate the arrays
	for (i = 0; i < TEST_SOURCES; i++) {
		if (posix_memalign(&buf, 64, TEST_LEN)) {
			printf("alloc error: Fail");
			return -1;
		}
		buffs[i] = buf;
	}

	for (i = 0; i < TEST_SOURCES; i++) {
		if (posix_memalign(&buf, 64, TEST_LEN)) {
			printf("alloc error: Fail");
			return -1;
		}
		temp_buffs[i] = buf;
	}

	// Test erasure code by encode and recovery

	encode_matrix = malloc(MMAX * KMAX);
	decode_matrix = malloc(MMAX * KMAX);
	invert_matrix = malloc(MMAX * KMAX);
	g_tbls = malloc(KMAX * TEST_SOURCES * 32);
	if (encode_matrix == NULL ||decode_matrix == NULL
		|| invert_matrix == NULL || g_tbls == NULL) {
		printf("Test failure! Error with malloc\n");
		return -1;
	};


	m = 9;
	k = 6;

	if (m > MMAX || k > KMAX)
		return -1; //exit

	//clock_t tic = clock();
	//clock_t start = clock();

	//printf("RANDOM DATA: \n");
	for (i = 0; i < k; i++){ //i is the dataword lengths
		for (j = 0; j < TEST_LEN; j++)
		{
			buffs[i][j] = rand();

			//printf("%u ", buffs[i][j]);
			//if ((j % TEST_LEN) == TEST_LEN-1)
			//printf("\n");
		}
	}

	//for (i = 0; i < k; i++){
	//	printf("%u ", *buffs[i]);
	//}


	//printf("\n");
	// Generate encode matrix encode_matrix
	// The matrix generated by gf_gen_rs_matrix
	gf_gen_cauchy1_matrix(encode_matrix, m, k);
	//output: encode_matrix

	/*printf("GENERATOR MATRIX: \n");
	for(i = 0 ; i < k*m ; i++)
	{
		printf("%u ",encode_matrix[i]);
		if((i % k) == k-1)
		printf("\n");
	}
	printf("\n");*/


	// Generate g_tbls from encode matrix
	ec_init_tables(k, m - k, &encode_matrix[k * k], g_tbls);
	// gf_vect_mul_init: Calculates const table gftbl in GF(2^8) from single input A
	//m - k = rows = the number of output vector, e.g., 9 - 5 = 4
	//&encode_matrix[k * k]: Pointer to sets of arrays of input coefficients used to encode or decode data
	//g_tbls: Pointer to start of space for concatenated output tables
 	//generated from input coefficients.  Must be of size 32*k*rows

	/*printf("Initial table: \n");
	for (i = 0; i < (32*(k*(m-k))); i++)
	{
		printf("%u ", g_tbls[i]);
		if ((i % 32) == 31)
		printf("\n");
	}*/

	// Perform matrix dot_prod for EC encoding
	// using g_tbls from encode matrix encode_matrix

	ec_encode_data(TEST_LEN, k, m - k, g_tbls, buffs, &buffs[k]);

	/* TEST_LEN: Length of each block of data (vector) of source data
	k: The number of vector source in G matrix
	m-k: The number of output vectors to concurrently encode/decode.
	gftbls: Pointer to array of input tables generated from coding
	coefficients in ec_init_tables(). Must be of size 32*k*rows
	buffs (data): Array of pointers to source input buffers
	&buffs[k]: Arrays of pointer to coded output buffers*/
	/*printf("\n");
	printf("ENCODED DATA: \n");
	for (i = 0; i < m; i++){
		for (j = 0; j < TEST_LEN; j++){
			printf("%u ", buffs[i][j]);
			if((j % TEST_LEN ) == TEST_LEN-1)
			printf("\n");
		}
	}
	*/


	// Choose random buffers to be in erasure
	memset(src_in_err, 0, TEST_SOURCES);
	//Fill whole src_in_err array with 0
	//TEST_SOURCE: Number of byte to be filled starting from src_in_err to be filled

	//Add errorsposix_memalign
	gen_err_list(src_err_list, src_in_err, &nerrs, &nsrcerrs, k, m);

	//BEGINNING OF TRACE REPAIR

	clock_t start = clock();
	precompute();

	unsigned char* rec;
	//unsigned char rec;
	unsigned char org[nerrs][TEST_LEN];
	int a, b;

	//original data
	/*
	for (a = 0; a < nerrs; a++){
		for (b = 0; b < TEST_LEN; b++){
			org[a][b] = buffs[src_err_list[a]][b];
		}
	} */

	for(int err = 0; err < nerrs; err++){
		j = src_err_list[err];
		rec = repair_trace_optimised(m, j, buffs);

		//compare recovered data with original data


		//this is check step for repair_trace_optimised
		for (b = 0; b < TEST_LEN; b++){
			if (rec[b] != buffs[j][b]) {
				printf("ERROR: TraceRepair doesn't match\n");
				exit(1);
			}
		}

		/*
		//This for repair_trace instead of repair_trace_optimised
		//unsigned char rec;
		for(int test_codeword = 0; test_codeword < TEST_LEN; test_codeword++){
			rec = repair_trace(9, 9, 9, j, buffs, test_codeword);
			//if (rec != buffs[j][test_codeword]){
			//	printf("ERROR: TraceRepair doesn't match\n");
			//	return -1;
			//}
		}*/
	}
	clock_t finish = clock();

	printf("-----------------------------------------\n");
	double trace_repair_time = ((double)(finish - start) / CLOCKS_PER_SEC)/3.0;
    printf("Elapsed (Repair_scheme): %f seconds\n", trace_repair_time);
	//printf("Run time (Rpair Trace): %d micro seconds\n",elapsed);
	printf("-----------------------------------------\n");
	//END OF TRACE REPAIR

	/*
	start = clock();
	unsigned char v;
	for (int iter = 0; iter < 300000; iter++){
		for (i = 0; i < 9; i++){
			for (j = 0; j < 9; j++){
				for (int u = 0; u < 72; u++){
					v = h_htbl[i][j][u];
				}
			}
		}
	}
	finish = clock();
	printf("Read table: %.10f seconds\n", 1.0*((double)(finish - start) / CLOCKS_PER_SEC));
	*/

	clock_t tic = clock();
	// Generate decode matrix from encode matrix
	re = gf_gen_decode_matrix(encode_matrix, decode_matrix,
				  invert_matrix, decode_index, src_err_list, src_in_err,
				  nerrs, nsrcerrs, k, m);
	if (re != 0) {
		printf("Fail to gf_gen_decode_matrix\n");
		return -1;
	}

	//printf("\n");
	/*
	//Print error list and non-erasure index
	printf(" - erase list = ");
	for (j = 0; j < nerrs; j++)
		printf(" %d", src_err_list[j]);
	printf("\n - Index = ");
	for (p = 0; p < k; p++){
		printf(" %d", decode_index[p]);
	}*/
	/*
	printf("\n \n");
	printf("nsrcerrs (number of original data error): %u \n", nsrcerrs);
	//mean the data are erasured that is in original data. Eg: lost data where<k.
	printf("nerrs (number of erased data): %u \n", nerrs);
	printf("\n");


	//Generate invert matrix
	printf("INVERT MATRIX: \n");
	for (i = 0; i < k*k; i++){
		printf(" %u", invert_matrix[i]);
		if ((i%k) == k-1)
		printf("\n");
	}
	printf("\n");*/


	//Generate decode matrix
	/*printf("DECODE MATRIX: \n");
	for (i = 0; i < nsrcerrs; i++){
		for(j = 0; j <k; j++ ){
			printf("%u ", decode_matrix[k*i + j]);
		}
		printf("\n");
	}
	for (p = nsrcerrs; p < nerrs; p++) { //nserrs: number of erased data unit , nerrs: Erased index
		for (i = 0; i < k; i++) {
			printf("%u ", decode_matrix[p * k + i]);
			if ((i%k) == k-1)
			printf("\n");
		}
	}
	printf("\n");*/


	// Pack recovery array as list of valid sources
	// Its order must be the same as the order
	// to generate matrix b in gf_gen_decode_matrix

	//printf("Recovery index: ");
	for (i = 0; i < k; i++) {
		recov[i] = buffs[decode_index[i]];
		//printf("%u ", *recov[i]);
	}

	//printf("\n \n");

	// Recover data
	ec_init_tables(k, nerrs, decode_matrix, g_tbls); //Generate or decode erasure codes on blocks of data
	/*printf("Initial table (recover part): \n");
	for (i = 0; i < 32*(k*(m-k)); i++)
	{
		printf("%u ", g_tbls[i]);
		if ((i % 32) == 31)
		printf("\n");
	}
	printf("\n");*/

	ec_encode_data(TEST_LEN, k, nerrs, g_tbls, recov, &temp_buffs[k]);// nerrs = rows, recov = data,
	//gettimeofday(&bye,NULL);
	clock_t toc = clock();
	double ISAL_repair_time = (double)(toc - tic) / CLOCKS_PER_SEC / 3.0;
	printf("Elapsed (ISA-L): %f seconds\n", ISAL_repair_time);
	//int elapsed1 = (((bye.tv_sec - hi.tv_sec) * 1000000) + (bye.tv_usec - hi.tv_usec)) / TEST_LEN;
	//printf("Run time (ISA-L): %d micro seconds\n",elapsed1);
	printf("-----------------------------------------\n");

	printf("Slow Factor: %.0f \n", trace_repair_time/ISAL_repair_time);
	/*for (i = 0; i < nerrs; i++){
		printf("RECOVER %u: \n", src_err_list[i]);
		dump(temp_buffs[k + i], TEST_LEN);
	}
	printf("\n");*/






	/*for (i = 0; i < nerrs; i++) {

		if (0 != memcmp(temp_buffs[k + i], buffs[src_err_list[i]], TEST_LEN)) { //compare size of temp_buff and buff
			printf("Fail error recovery (%d, %d, %d)\n", m, k, nerrs);
			printf(" - erase list = ");
			for (j = 0; j < nerrs; j++)
				printf(" %d", src_err_list[j]);
			printf(" - Index = ");
			for (p = 0; p < k; p++)
				printf(" %d", decode_index[p]);
			printf("\nencode_matrix:\n");
			dump_u8xu8((u8 *) encode_matrix, m, k);
			printf("inv b:\n");
			dump_u8xu8((u8 *) invert_matrix, k, k);
			printf("\ndecode_matrix:\n");
			dump_u8xu8((u8 *) decode_matrix, m, k);
			printf("recov %d:", src_err_list[i]);
			dump(temp_buffs[k + i], k*m);
			printf("orig   :");
			dump(buffs[src_err_list[i]], k*m);
			return -1;
		}
	}*/

	//printf("\nEnter j: ");
	//scanf("%d", &j);

	/*m = 9;
	k = 5;
	if (m > MMAX || k > KMAX)
		return -1;

	// Make random data
	for (i = 0; i < k; i++)
		for (j = 0; j < TEST_LEN; j++)
			buffs[i][j] = rand();

	// The matrix generated by gf_gen_cauchy1_matrix
	// is always invertable.
	gf_gen_cauchy1_matrix(encode_matrix, m, k);

	// Generate g_tbls from encode matrix encode_matrix
	ec_init_tables(k, m - k, &encode_matrix[k * k], g_tbls);

	// Perform matrix dot_prod for EC encoding
	// using g_tbls from encode matrix encode_matrix
	ec_encode_data(TEST_LEN, k, m - k, g_tbls, buffs, &buffs[k]);  //&buff[k] = encode_matrix

	// Choose random buffers to be in erasure
	memset(src_in_err, 0, TEST_SOURCES);
	gen_err_list(src_err_list, src_in_err, &nerrs, &nsrcerrs, k, m);

	// Generate decode matrix
	re = gf_gen_decode_matrix(encode_matrix, decode_matrix,
				  invert_matrix, decode_index, src_err_list, src_in_err,
				  nerrs, nsrcerrs, k, m);
	if (re != 0) {
		printf("Fail to gf_gen_decode_matrix\n");
		return -1;
	}
	// Pack recovery array as list of valid sources
	// Its order must be the same as the order
	// to generate matrix b in gf_gen_decode_matrix
	for (i = 0; i < k; i++) {
		recov[i] = buffs[decode_index[i]];
	}

	// Recover data
	ec_init_tables(k, nerrs, decode_matrix, g_tbls);
	ec_encode_data(TEST_LEN, k, nerrs, g_tbls, recov, &temp_buffs[k]);
	for (i = 0; i < nerrs; i++) {

		if (0 != memcmp(temp_buffs[k + i], buffs[src_err_list[i]], TEST_LEN)) {
			printf("Fail error recovery (%d, %d, %d)\n", m, k, nerrs);
			printf(" - erase list = ");
			for (j = 0; j < nerrs; j++)
				printf(" %d", src_err_list[j]);
			printf(" - Index = ");
			for (p = 0; p < k; p++)
				printf(" %d", decode_index[p]);
			printf("\nencode_matrix:\n");
			dump_u8xu8((u8 *) encode_matrix, m, k);
			printf("inv b:\n");
			dump_u8xu8((u8 *) invert_matrix, k, k);
			printf("\ndecode_matrix:\n");
			dump_u8xu8((u8 *) decode_matrix, m, k);
			printf("recov %d:", src_err_list[i]);
			dump(temp_buffs[k + i], 25);
			printf("orig   :");
			dump(buffs[src_err_list[i]], 25);
			return -1;
		}
	}

	// Do more random tests
	for (rtest = 0; rtest < RANDOMS; rtest++) {
		while ((m = (rand() % MMAX)) < 2) ;
		while ((k = (rand() % KMAX)) >= m || k < 1) ;

		if (m > MMAX || k > KMAX)
			continue;

		// Make random data
		for (i = 0; i < k; i++)
			for (j = 0; j < TEST_LEN; j++)
				buffs[i][j] = rand();

		// The matrix generated by gf_gen_cauchy1_matrix
		// is always invertable.
		gf_gen_cauchy1_matrix(encode_matrix, m, k);

		// Make parity vects
		// Generate g_tbls from encode matrix a
		ec_init_tables(k, m - k, &encode_matrix[k * k], g_tbls);
		// Perform matrix dot_prod for EC encoding
		// using g_tbls from encode matrix a
		ec_encode_data(TEST_LEN, k, m - k, g_tbls, buffs, &buffs[k]);

		// Random errors
		memset(src_in_err, 0, TEST_SOURCES);
		gen_err_list(src_err_list, src_in_err, &nerrs, &nsrcerrs, k, m);

		// Generate decode matrix
		re = gf_gen_decode_matrix(encode_matrix, decode_matrix,
					  invert_matrix, decode_index, src_err_list,
					  src_in_err, nerrs, nsrcerrs, k, m);
		if (re != 0) {
			printf("Fail to gf_gen_decode_matrix\n");
			return -1;
		}
		// Pack recovery array as list of valid sources
		// Its order must be the same as the order
		// to generate matrix b in gf_gen_decode_matrix
		for (i = 0; i < k; i++) {
			recov[i] = buffs[decode_index[i]];
		}

		// Recover data
		ec_init_tables(k, nerrs, decode_matrix, g_tbls);
		ec_encode_data(TEST_LEN, k, nerrs, g_tbls, recov, &temp_buffs[k]);

		for (i = 0; i < nerrs; i++) {

			if (0 != memcmp(temp_buffs[k + i], buffs[src_err_list[i]], TEST_LEN)) {
				printf("Fail error recovery (%d, %d, %d) - ", m, k, nerrs);
				printf(" - erase list = ");
				for (j = 0; j < nerrs; j++)
					printf(" %d", src_err_list[j]);
				printf(" - Index = ");
				for (p = 0; p < k; p++)
					printf(" %d", decode_index[p]);
				printf("\nencode_matrix:\n");
				dump_u8xu8((u8 *) encode_matrix, m, k);
				printf("inv b:\n");
				dump_u8xu8((u8 *) invert_matrix, k, k);
				printf("\ndecode_matrix:\n");
				dump_u8xu8((u8 *) decode_matrix, m, k);
				printf("orig data:\n");
				dump_matrix(buffs, m, 25);
				printf("orig   :");
				dump(buffs[src_err_list[i]], 25);
				printf("recov %d:", src_err_list[i]);
				dump(temp_buffs[k + i], 25);
				return -1;
			}
		}
		putchar('.');
	}

	// Run tests at end of buffer for Electric Fence
	k = 16;
	align = (LEN_ALIGN_CHK_B != 0) ? 1 : 16;
	if (k > KMAX)
		return -1;

	for (rows = 1; rows <= 16; rows++) {
		m = k + rows;
		if (m > MMAX)
			return -1;

		// Make random data
		for (i = 0; i < k; i++)
			for (j = 0; j < TEST_LEN; j++)
				buffs[i][j] = rand();

		for (size = EFENCE_TEST_MIN_SIZE; size <= EFENCE_TEST_MAX_SIZE; size += align) {
			for (i = 0; i < m; i++) {	// Line up TEST_SIZE from end
				efence_buffs[i] = buffs[i] + TEST_LEN - size;
			}

			// The matrix generated by gf_gen_cauchy1_matrix
			// is always invertable.
			gf_gen_cauchy1_matrix(encode_matrix, m, k);

			// Make parity vects
			// Generate g_tbls from encode matrix a
			ec_init_tables(k, m - k, &encode_matrix[k * k], g_tbls);
			// Perform matrix dot_prod for EC encoding
			// using g_tbls from encode matrix a
			ec_encode_data(size, k, m - k, g_tbls, efence_buffs, &efence_buffs[k]);

			// Random errors
			memset(src_in_err, 0, TEST_SOURCES);
			gen_err_list(src_err_list, src_in_err, &nerrs, &nsrcerrs, k, m);

			// Generate decode matrix
			re = gf_gen_decode_matrix(encode_matrix, decode_matrix,
						  invert_matrix, decode_index, src_err_list,
						  src_in_err, nerrs, nsrcerrs, k, m);
			if (re != 0) {
				printf("Fail to gf_gen_decode_matrix\n");
				return -1;
			}
			// Pack recovery array as list of valid sources
			// Its order must be the same as the order
			// to generate matrix b in gf_gen_decode_matrix
			for (i = 0; i < k; i++) {
				recov[i] = efence_buffs[decode_index[i]];
			}

			// Recover data
			ec_init_tables(k, nerrs, decode_matrix, g_tbls);
			ec_encode_data(size, k, nerrs, g_tbls, recov, &temp_buffs[k]);

			for (i = 0; i < nerrs; i++) {

				if (0 !=
				    memcmp(temp_buffs[k + i], efence_buffs[src_err_list[i]],
					   size)) {
					printf("Efence: Fail error recovery (%d, %d, %d)\n", m,
					       k, nerrs);

					printf("size = %d\n", size);

					printf("Test erase list = ");
					for (j = 0; j < nerrs; j++)
						printf(" %d", src_err_list[j]);
					printf(" - Index = ");
					for (p = 0; p < k; p++)
						printf(" %d", decode_index[p]);
					printf("\nencode_matrix:\n");
					dump_u8xu8((u8 *) encode_matrix, m, k);
					printf("inv b:\n");
					dump_u8xu8((u8 *) invert_matrix, k, k);
					printf("\ndecode_matrix:\n");
					dump_u8xu8((u8 *) decode_matrix, m, k);

					printf("recov %d:", src_err_list[i]);
					dump(temp_buffs[k + i], align);
					printf("orig   :");
					dump(efence_buffs[src_err_list[i]], align);
					return -1;
				}
			}
		}

	}

	// Test rand ptr alignment if available

	for (rtest = 0; rtest < RANDOMS; rtest++) {
		while ((m = (rand() % MMAX)) < 2) ;
		while ((k = (rand() % KMAX)) >= m || k < 1) ;

		if (m > MMAX || k > KMAX)
			continue;

		size = (TEST_LEN - PTR_ALIGN_CHK_B) & ~15;

		offset = (PTR_ALIGN_CHK_B != 0) ? 1 : PTR_ALIGN_CHK_B;
		// Add random offsets
		for (i = 0; i < m; i++) {
			memset(buffs[i], 0, TEST_LEN);	// zero pad to check write-over
			memset(temp_buffs[i], 0, TEST_LEN);	// zero pad to check write-over
			ubuffs[i] = buffs[i] + (rand() & (PTR_ALIGN_CHK_B - offset));
			temp_ubuffs[i] = temp_buffs[i] + (rand() & (PTR_ALIGN_CHK_B - offset));
		}

		for (i = 0; i < k; i++)
			for (j = 0; j < size; j++)
				ubuffs[i][j] = rand();

		// The matrix generated by gf_gen_cauchy1_matrix
		// is always invertable.
		gf_gen_cauchy1_matrix(encode_matrix, m, k);

		// Make parity vects
		// Generate g_tbls from encode matrix a
		ec_init_tables(k, m - k, &encode_matrix[k * k], g_tbls);
		// Perform matrix dot_prod for EC encoding
		// using g_tbls from encode matrix a
		ec_encode_data(size, k, m - k, g_tbls, ubuffs, &ubuffs[k]);

		// Random errors
		memset(src_in_err, 0, TEST_SOURCES);
		gen_err_list(src_err_list, src_in_err, &nerrs, &nsrcerrs, k, m);

		// Generate decode matrix
		re = gf_gen_decode_matrix(encode_matrix, decode_matrix,
					  invert_matrix, decode_index, src_err_list,
					  src_in_err, nerrs, nsrcerrs, k, m);
		if (re != 0) {
			printf("Fail to gf_gen_decode_matrix\n");
			return -1;
		}
		// Pack recovery array as list of valid sources
		// Its order must be the same as the order
		// to generate matrix b in gf_gen_decode_matrix
		for (i = 0; i < k; i++) {
			recov[i] = ubuffs[decode_index[i]];
		}

		// Recover data
		ec_init_tables(k, nerrs, decode_matrix, g_tbls);
		ec_encode_data(size, k, nerrs, g_tbls, recov, &temp_ubuffs[k]);

		for (i = 0; i < nerrs; i++) {

			if (0 != memcmp(temp_ubuffs[k + i], ubuffs[src_err_list[i]], size)) {
				printf("Fail error recovery (%d, %d, %d) - ", m, k, nerrs);
				printf(" - erase list = ");
				for (j = 0; j < nerrs; j++)
					printf(" %d", src_err_list[j]);
				printf(" - Index = ");
				for (p = 0; p < k; p++)
					printf(" %d", decode_index[p]);
				printf("\nencode_matrix:\n");
				dump_u8xu8((unsigned char *)encode_matrix, m, k);
				printf("inv b:\n");
				dump_u8xu8((unsigned char *)invert_matrix, k, k);
				printf("\ndecode_matrix:\n");
				dump_u8xu8((unsigned char *)decode_matrix, m, k);
				printf("orig data:\n");
				dump_matrix(ubuffs, m, 25);
				printf("orig   :");
				dump(ubuffs[src_err_list[i]], 25);
				printf("recov %d:", src_err_list[i]);
				dump(temp_ubuffs[k + i], 25);
				return -1;
			}
		}

		// Confirm that padding around dests is unchanged
		memset(temp_buffs[0], 0, PTR_ALIGN_CHK_B);	// Make reference zero buff

		for (i = 0; i < m; i++) {

			offset = ubuffs[i] - buffs[i];

			if (memcmp(buffs[i], temp_buffs[0], offset)) {
				printf("Fail rand ualign encode pad start\n");
				return -1;
			}
			if (memcmp
			    (buffs[i] + offset + size, temp_buffs[0],
			     PTR_ALIGN_CHK_B - offset)) {
				printf("Fail rand ualign encode pad end\n");
				return -1;
			}
		}

		for (i = 0; i < nerrs; i++) {

			offset = temp_ubuffs[k + i] - temp_buffs[k + i];
			if (memcmp(temp_buffs[k + i], temp_buffs[0], offset)) {
				printf("Fail rand ualign decode pad start\n");
				return -1;
			}
			if (memcmp
			    (temp_buffs[k + i] + offset + size, temp_buffs[0],
			     PTR_ALIGN_CHK_B - offset)) {
				printf("Fail rand ualign decode pad end\n");
				return -1;
			}
		}

		putchar('.');
	}

	// Test size alignment

	align = (LEN_ALIGN_CHK_B != 0) ? 13 : 16;

	for (size = TEST_LEN; size > 0; size -= align) {
		while ((m = (rand() % MMAX)) < 2) ;
		while ((k = (rand() % KMAX)) >= m || k < 1) ;

		if (m > MMAX || k > KMAX)
			continue;

		for (i = 0; i < k; i++)
			for (j = 0; j < size; j++)
				buffs[i][j] = rand();

		// The matrix generated by gf_gen_cauchy1_matrix
		// is always invertable.
		gf_gen_cauchy1_matrix(encode_matrix, m, k);

		// Make parity vects
		// Generate g_tbls from encode matrix a
		ec_init_tables(k, m - k, &encode_matrix[k * k], g_tbls);
		// Perform matrix dot_prod for EC encoding
		// using g_tbls from encode matrix a
		ec_encode_data(size, k, m - k, g_tbls, buffs, &buffs[k]);

		// Random errors
		memset(src_in_err, 0, TEST_SOURCES);
		gen_err_list(src_err_list, src_in_err, &nerrs, &nsrcerrs, k, m);
		// Generate decode matrix
		re = gf_gen_decode_matrix(encode_matrix, decode_matrix,
					  invert_matrix, decode_index, src_err_list,
					  src_in_err, nerrs, nsrcerrs, k, m);
		if (re != 0) {
			printf("Fail to gf_gen_decode_matrix\n");
			return -1;
		}
		// Pack recovery array as list of valid sources
		// Its order must be the same as the order
		// to generate matrix b in gf_gen_decode_matrix
		for (i = 0; i < k; i++) {
			recov[i] = buffs[decode_index[i]];
		}

		// Recover data
		ec_init_tables(k, nerrs, decode_matrix, g_tbls);
		ec_encode_data(size, k, nerrs, g_tbls, recov, &temp_buffs[k]);

		for (i = 0; i < nerrs; i++) {

			if (0 != memcmp(temp_buffs[k + i], buffs[src_err_list[i]], size)) {
				printf("Fail error recovery (%d, %d, %d) - ", m, k, nerrs);
				printf(" - erase list = ");
				for (j = 0; j < nerrs; j++)
					printf(" %d", src_err_list[j]);
				printf(" - Index = ");
				for (p = 0; p < k; p++)
					printf(" %d", decode_index[p]);
				printf("\nencode_matrix:\n");
				dump_u8xu8((unsigned char *)encode_matrix, m, k);
				printf("inv b:\n");
				dump_u8xu8((unsigned char *)invert_matrix, k, k);
				printf("\ndecode_matrix:\n");
				dump_u8xu8((unsigned char *)decode_matrix, m, k);
				printf("orig data:\n");
				dump_matrix(buffs, m, 25);
				printf("orig   :");
				dump(buffs[src_err_list[i]], 25);
				printf("recov %d:", src_err_list[i]);
				dump(temp_buffs[k + i], 25);
				return -1;
			}
		}
	}

	printf("done EC tests: Pass\n");*/
	return 0;
}
