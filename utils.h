// some utility functions

// linear interpolation of the function f(xs) = ys
// at point xi which must be min(xs) <= xi <= max(xs)
double linear_interp(const double xi, const double* xs, const double* ys);
