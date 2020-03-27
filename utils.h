// some utility functions

// linear interpolation of the function f(xs) = ys
// at point xi which must be min(xs) <= xi <= max(xs)
double linear_interp(const double xi, const double* xs, const double* ys) {
	double m,b,ya,yb,xa,xb,yc;
	int i=0;

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