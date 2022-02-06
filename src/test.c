#include <stdio.h>
#include <stdlib.h>
#include <string.h>		// for memset, memcmp
#include "erasure_code.h"

int main()
{
    int m;
    int k;
    unsigned char *encode_matrix;
    unsigned char **s;
    m = 9; 
    k = 6;
    gf_gen_cauchy1_matrix(encode_matrix, m, k);
        int i, j;
	    for (i = 0; i < k; i++) 
        {
		    for (j = 0; j < m; j++) 
            {
			    printf(" %2x", s[m][k]);
		    }
		    printf("\n");
	    }
	printf("\n");
    return 0;
}  