# Trace Repair for Reed-Solomon Codes
Implement the trace repair method to recover a single erasure for Reed-Solomon codes on top of ISA-L  
Contact sonhoang.dau@rmit.edu.au for any feedback or questions  
  
The main trace repair method is implemented in  
        unsigned char* repair_trace_optimised(int n, int j, unsigned char **buffs, double* repair_time_ptr)  
in which  
-n: code length, e.g., n = 9, then the codeword c = (c0,c1,...,c8)  
-j: an integer between 0 and n-1 representing the index to be recovered (cj is lost)  
-buffs: a 2-D array of dimension nxTEST_LEN, where TEST_LEN = 10,000,000 is the number of test codewords  
        buffs is generated randomly in main()  
-repair_time_ptr: record the total running time  
   
# Parameters Setup  
-Set TEST_LEN = 10,000,000 for example  
-Go to main(), set m = 9, k = 6 for a RS(9,6), note that the code use m in main(), not n  
-Go to the preamble, uncomment the corresponding header file to access the required tables H, R, D, e.g.,   
        #include "repair_scheme_tbl_9_6_ISIT.h"  
  
# Compilation & Run  
-Normal compilation and run on IDE  
-Linux:  
        gcc -std=gnu99 erasure_code_test_optimised_ISIT.c -o erasure_code_test_optimised_ISIT  
        ./erasure_code_test_optimised_ISIT  
