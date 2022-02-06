#include <limits.h>
#include <string.h>		// for memset
#include "erasure_code.h"
int main()
{
    
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
}
