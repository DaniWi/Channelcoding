#include <Rcpp.h>
using namespace Rcpp;

int getBinaryLength(int number) {
	
	int bits = 0;
	
	while (number >>= 1) {
	    bits++;
	}
	
	return bits;
}

int polynom_division_remainder_mod2(int numerator, int denominator) {
	// int quotient = 0;
	
	int numeratorLength = getBinaryLength(numerator);
	int denominatorLength = getBinaryLength(denominator);
	
	int iterations = numeratorLength - denominatorLength + 1;
	
	for (int i = 0; i < iterations; i++) {
		
		if ((numerator & (int)pow(2,numeratorLength-i)) != 0) {
			
			// the following line computes the quotient but is not of interes in this application
			// quotient += pow(2,numeratorLength - denominatorLength - i);
			
			numerator ^= denominator << (numeratorLength - denominatorLength - i);
		}
	}
	
	return numerator;	// matches the remainder
}

int gcd_polynomial(int a, int b) {
	if (b > a) {
		int c = a;
		a = b;
		b = c;
	}

	int r = polynom_division_remainder_mod2(a,b);
	while (r != 0) {
		a = b;
		b = r;
		r = polynom_division_remainder_mod2(a,b);
	}
	
	return b;
	
	/*
	while (b > 0) {
		int t = a;
		a = b;
		b = polynom_division_remainder_mod2(t,b);
	}
	
	return a;
	*/
}

// [[Rcpp::export]]
int gcd_polynomial(IntegerVector x) {
	
	int gcd = 0;
	
	int vector_size = x.size();
	
	if (vector_size == 1) {
		gcd = x[1];
	} else if (vector_size == 2) {
		gcd = gcd_polynomial(x[0], x[1]);
	} else {
		gcd = gcd_polynomial(x[0], x[1]);
	    for (int i = 2; i < vector_size; i++) {
	    	gcd = gcd_polynomial(gcd, x[i]);
	    	if (gcd == 1) { break; }
		}
	}
	
	return gcd;
}
