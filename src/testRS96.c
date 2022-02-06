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

	return gff_base_f16[(i = gflog_base_f16[a] + gflog_base_f16[b]) > 14 ? i - 15 : i];
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
	/*unsigned char a[k*m];
	for (i = 0; i < k*m; i++)
		a[i] = 0;

	for (i = 0; i < k; i++)
		a[k * i + i] = 1;

	for (i = k; i < m; i++) 
	{
		p = 1;
	
		for (j = 0; j < k; j++) 
		{
			a[k * i + j] = p;
			p = gf_mul(p, gen);
		}
		
		gen = gf_mul(gen, 2);
	
	}*/

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

int main(int argc, char * argv[])
{
    int i, j;
	int m = 9;
	int k = 6;

	unsigned char *invert_matrix;
	unsigned char *encode_matrix;
	int n = 6;

	encode_matrix = malloc(127 * 127);
	invert_matrix = malloc(127 * 127);

	gf_gen_rs_matrix(encode_matrix, m, k);

	gf_invert_matrix(&encode_matrix[k*m], invert_matrix, n);

	printf("Encode matrix: \n");
	for (i = 0; i < k*m; i++){
		printf("%u ", encode_matrix[i]);
		if ((i%k) == k-1)
		printf("\n");

	}
	printf("\n");

	printf("Inverse matrix: \n");
	for (i = 0; i < k*k; i++){
		printf(" %u", invert_matrix[i]);
		if ((i%k) == k-1)
		printf("\n");
	}
	printf("\n");
	
	

	/*for(i = 0 ; i < k*m ; i++)
	{
		printf("%u ",a[i]);
		if((i % k) == k-1)
		printf("\n");
	}*/

	

	return 0;
}



