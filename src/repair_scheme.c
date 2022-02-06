#include <stdio.h>
#include <stdlib.h>
#include <string.h>	
#include<math.h>
#include <limits.h>	
#include "erasure_code.h"
#include "types.h"
#include "ec_base.h"
#include "repair_scheme_tbl.h"

#define KMAX 6
#define MMAX 9
#define TEST_SOURCES 127
#define TEST_LEN 1
#define l 8
#define BIN_LEN 8


#define swap(head, tail) \
		unsigned char temp; \
		temp = head; \
		head = tail; \
		tail = temp;	


typedef unsigned char u8;
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

void gf_gen_cauchy_matrix(unsigned char *a, int m, int k)
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
{
	int i, j, k;
	unsigned char s;
	
	for (k = 0; l < dests; k++) { // l = 1
		for (i = 0; i < len; i++) { // i = 0
			s = 0;
			for (j = 0; j < srcs; j++){
				s^= gf_mul(src[j][i], v[j * 32 + l * srcs * 32 + 1]);
			}
			dest[l][i] = s;
			printf("s: %u \n", s);
		}
	}
}


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

/***************************************************************/
// Convert decimal to binary 
int binRep(unsigned char N) 
{ 
    // To store the binary number 
    int B_Number = 0; 
    int c;
    int count = 0; 

    while (N != 0) { 
        int rem = N % 2; 
        c = pow(10, count); 

        B_Number += rem * c; 

        N /= 2; 
        count++; 
    }
    return B_Number;
}


//Decalre traceSum variable as 1 bit
struct{
	unsigned char traceSum:1;
	unsigned char rh:1;

}bit_field;

//Input m and return t-bit (array). For ColumnTrc calculation
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


unsigned char binSum( int x, int y)
{ 
    int i;
	unsigned char r1[8], r2[8], r[32];
	unsigned char *p1, *p2;
	
	p1 = malloc(8);
	p1 =  binGap(8, x);

	int tail = 7;

    for ( i = 0; i < 8; i++)
    {
        r1[i] = p1[i];
    }

	for ( i = 0; i < 4; i++)
	{
		swap(r1[i], r1[tail]);
		tail--;
	}

    p2 = malloc(8);
	p2 = binGap(8, y);

    for ( i = 0; i < 8; i++)
    {
        r2[i] = p2[i];
		//printf("%u", r1[i]);
    }

	bit_field.traceSum = 0;
	for(i = 0; i < l; i++){
		bit_field.traceSum ^= (r1[i]*r2[i]);
	}
	
	free(p1);
	free(p2);
	
    return bit_field.traceSum;
}

      

//MAIN FUNCTION: Helper Node Banwidth calculation and Generate Traces from Helper Node to send to Recover Node
unsigned char repair_trace(int block, int row, int column, int j, unsigned char **buffs, 
					int htbl[][9][9],int rtbl[][9][9],int dtbl[][8]){

	int a, b; // for FOR loop 
    unsigned int bw;
	int s, vi; 
	int i;
	unsigned char bi[MMAX];
	unsigned char *vec;
	unsigned char RepairTr[MMAX][l];
	unsigned char ColumnTr[MMAX][l];
	unsigned char TargetTr[MMAX];
	unsigned char rev, trace;

    //H stands for Helper node, R stands for recover node
	vec = malloc(MMAX * MMAX);
    bw = 0;

    //Calculate bandwidth using H table
	printf("Bandwidth: ");
    for (i = 0; i < block; i++){
		if (i != j){
			bw = htbl[i][j][0];
		}
		else{
			continue;
		}
		bi[i] = bw;
		printf("[%d]: %u  ",i, bi[i]);
    }
	
	
	//print Cj
	for(i = 0; i <block; i++){
		printf("\nC_j[%d]: %u",i, *buffs[i]);
	}
	printf("\n");

	//Calculate traces to send to Node j (H table)
	printf("\nRepair Traces: \n");
    for (i = 0; i < block ; i++){ 
		if( i != j){
        	for (a = 0; a < bi[i]; a++){
				RepairTr[i][a] = binSum(htbl[i][j][a+1], (int)*buffs[i]);
				//RepairTr[i][a] = trace;
				printf(" %u ", RepairTr[i][a]);
			}
		}
		else
		{
			continue;
		}
		printf("\n");
    }
	printf("\n");


	//Recover node generates Column traces using R table
	printf("Column Tr: \n");
    for (i = 0; i < block; i++){
		if (i != j){
			unsigned char* p = RepairTr[i];
        	for (s = 0; s < l; s++){
				vi = rtbl[i][j][s]; 
				//printf(" vi : %u ", vi);
				vec = binGap(htbl[i][j][0], vi);

				//Print test BinRep 
				//printf(" Vec [%d] * P[a]: ", s);

				unsigned char result;
				result = 0;
				for(a = 0; a < bi[i]; a++){
					result ^= (vec[a]*p[a]);
					//printf(" (%u * %u) ^ ", vec[a], p[a]); // p[a] is the basis, vec is the coef vector
				}
				
				//printf("\n");
				ColumnTr[i][s] = result;
				printf(" %u ", ColumnTr[i][s]);
				//printf("\n");
			
				
			}
			printf("\n");
    	}
		else{
			continue;
		}
	}

	//Construct t Target Traces
	printf("\nTarget tr: \n");
	for (s = 0; s < l; s++){
		bit_field.rh = 0;//0*Z(q)
		for (i = 0; i < block; i++){
			if (i != j){
				bit_field.rh ^= ColumnTr[i][s];
				//printf(" %u ", ColumnTr[i][s]);
			}
			else{
				continue;
			}
		}
		//printf("\n");
		TargetTr[s] = bit_field.rh; 
		printf(" %u |",TargetTr[s]);
	} 

	printf("\n\n");
	//Assume that we obtained D table assigned as dtbl
	rev = 0; //0*Z(q)
	for(s = 0; s < l; s++){
		if(TargetTr[s] != 0){
			rev ^=dtbl[j][s];
			printf(" dtbl[%d]: %d ",s,dtbl[j][s]);
			
		}	
		else{
			continue;
		}
		printf("\n");
		
	}
	
	if(rev==(*buffs[j]))
		printf("\nTRUE");
	else
		printf("\nFALSE");
	
	free(vec);
	return rev;
}


int main(int argc, char *argv[])
{
	int i, j, m, k;
	int a, b; //FOR loop condition
	void *buf;
	unsigned char *buffs[TEST_SOURCES];
	unsigned char *encode_matrix, *g_tbls;
	int block, row, column;
	u8 *ubuffs[TEST_SOURCES];

    // Allocate the arrays
	for (a = 0; a < TEST_SOURCES; a++) {
		if (posix_memalign(&buf, 64, TEST_LEN)) {
			printf("alloc error: Fail");
			return -1; 
		}
		buffs[a] = buf;
	}


	encode_matrix = malloc(MMAX * KMAX);
	g_tbls = malloc(KMAX * TEST_SOURCES * 32);
	if (encode_matrix == NULL || g_tbls == NULL) {
		printf("Test failure! Error with malloc\n");
		return -1;
	};

    //STEP 1: Create codeword using generator matrix of RS(9, 6)
	// Pick a first test of RS(9, 6)
	m = 9;
	k = 6;
	if (m > MMAX || k > KMAX)
		return -1; //exi

	// Make random data	
	printf("RANDOM DATA: \n");
	for (a = 0; a < k; a++){ //i is the dataword lengths
		for (b = 0; b < TEST_LEN; b++) 
		{		
			scanf("%d", &buffs[a][b]);
			printf("%u ", buffs[a][b]);
			if ((b % TEST_LEN) == TEST_LEN-1)
			printf("\n");
		}
	}
	printf("\n");


	gf_gen_cauchy_matrix(encode_matrix, m, k); 
	printf("GENERATOR MATRIX: \n");
	for(a = 0 ; a < k*m ; a++)
	{
		printf("%u ",encode_matrix[a]);
		if((a % k) == k-1)
		printf("\n");
	}
	printf("\n");


	// Generate g_tbls from encode matrix
	ec_init_tables(k, m - k, &encode_matrix[k * k], g_tbls); 

    //encode data using init table
	ec_encode_data(TEST_LEN, k, m - k, g_tbls, buffs, &buffs[k]);
	printf("\n");
	printf("ENCODED DATA: \n");
	for (a = 0; a < m; a++){
		for (b = 0; b < TEST_LEN; b++){
			printf("%u ", buffs[a][b]);
			if((j % TEST_LEN ) == TEST_LEN-1)
			printf("\n");
		}
	}

	
	unsigned char buffs[9] = {221, 166, 119, 144, 48, 153, 241, 156, 5};
	//unsigned char buffs[9] = {41, 35, 190, 132, 225, 108, 56, 111, 68};

	//read file into array
    /*unsigned char buffs[9];
	FILE *myFile;
    myFile = fopen("RandomCodeword96.txt", "r");

    if (myFile == NULL){
        printf("Error Reading File\n");
        exit (0);
    }

    for (a = 0; a < 9; a++)
    {
        fscanf(myFile, "%d,", &buffs[a]);
    }
	fclose(myFile);*/

	//STEP 2: Test repair trace and bandwidth calculation
	printf("\nEnter j: ");
	scanf("%d", &j);
	
	printf("\nNode n0. %d is erased!\n\n", j);

	unsigned char rec;
	rec = repair_trace(9, 9, 9, j, buffs, h_htbl, h_rtbl, h_dtbl);
	printf("\nrecovery variable: %u", rec);

    return 0;
}


