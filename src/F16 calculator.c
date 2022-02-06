#include <stdio.h>
#include <stdlib.h>

unsigned long long int f2polydeg(unsigned long long int r)
{
	int d = 0;
	if (r == 0) {
		return 0;

	}
	while (r >>= 1)
		d++;
	return d;
}

unsigned long long int ff2mul(unsigned long long int a, unsigned long long int b, unsigned long long int r)
{
	int degb = f2polydeg(b);
	int degr = f2polydeg(r);
	unsigned long long int bit_at_degr = 1 << degr;
	unsigned long long int prod = 0;
	unsigned long long int temp;
	int i, power;

	for (i = 0; i <= degb; i++) {
		// Test if this bit position is set.
		if (!((b >> i) & 1))
			continue;

		// Multiply by the power of x on the current term of b, reducing mod r en route.
		temp = a;
		for (power = 1; power <= i; power++) {
			temp <<= 1;
			if (temp & bit_at_degr)
				temp ^= r;
		}
		// Add in this partial product.
		prod ^= temp;
	}

	return prod;
}

// Repeated-squaring algorithm.  This handles positive powers only.
unsigned int ff2power(unsigned int a, unsigned int power, unsigned int  r)
{
	unsigned a2 = a;
	unsigned out = 1;

	while (power != 0) {
		if (power & 1) // Test current bit
			out = ff2mul(out, a2, r);
		power = power >> 1; // Shift
		a2 = ff2mul(a2, a2, r);
	}
	return out;
}

unsigned long long int ff2recip(unsigned long long int b, unsigned long long int r)
{
	int n = f2polydeg(r);
	unsigned pn2 = (1 << n) - 2;
	return ff2power(b, pn2, r);
}

unsigned long long int ff2div(unsigned long long int a, unsigned long long int b, unsigned long long int r)
{
	unsigned long long int binv;
	binv = ff2recip(b, r);
	return ff2mul(a, binv, r);
}

int main(unsigned int a, unsigned int power, unsigned int r){
	a = 2;
	r = 285;
	
	printf("enter power: ");
	scanf("%u", &power);
	printf("%u ", ff2power(a, power, r));

	return 0;
}