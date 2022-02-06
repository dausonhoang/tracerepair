#include<stdlib.h>
#include<stdio.h>

int mul(int n, int m)
{
    int ans = 0, count = 0;
    while (m)
    {
        // check for set bit and left 
        // shift n, count times
        if (m % 2 == 1)              
            ans += n << count;
 
        // increment of place value (count)
        count++;
        m /= 2;
    }
    return ans;
}

int main(){

    int c; 
    c = mul(1, 1);
    printf("%d", c);
}