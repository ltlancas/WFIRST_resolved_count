// some utility functions
#include "nr3.h"


// linear interpolation of the function f(xs) = ys
// at point xi which must be min(xs) <= xi <= max(xs)
Doub linear_interp(const Doub xi, const Doub* xs, const Doub* ys) {
	Doub m,b,ya,yb,xa,xb,yc;
	Int i=0;

	// get to the right spot in the array
	while(xi<xs[i]) i++;
	i--;

	// linear interpolation 
	xa = xs[i];
	xb = xs[i+1];
	ya = ys[i];
	yb = ys[i+1];
	m = (yb-ya)/(xb-xa);
	b = yb - m*xb;
	yc = m*xi + b;

	return yc;
}